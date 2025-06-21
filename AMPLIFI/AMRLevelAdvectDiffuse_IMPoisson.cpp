//
//  AMRLevelAdvectDiffuse_IMPoisson.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 6/20/25.
//

#include <iostream>
#include <iomanip>

#include "parstream.H"

#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "computeNorm.H"
#include "OldTimer.H"

#include "AMRLevelAdvectDiffuse.H"
#include "AdvectPhysicsF_F.H"
#include "AdvectTestIBC.H"
#include "GodunovUtilitiesF_F.H"
#include "ParmParse.H"
#include "CellToEdge.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "EdgeToCell.H"
#include "AdvectPhysicsF_F.H"
#include "AdvectDiffuseUtils.H"
#include "AMRPoissonOp.H"
#include "VCAMRPoissonOp2.H"
#include "GodunovPhysicsF_F.H"
#include "NamespaceHeader.H"

#include "LevelTGAF_F.H"

#include "additionalBCFunc.hpp"
#include "photoionization.hpp"
#include "physicalConstants.h"
#include "Gradient.H"
#include "dataFileIFReduced.hpp"


/*******/
void
AMRLevelAdvectDiffuse::
poissonSolveImplicitComposite1() {
  if (m_hasFiner) {
    if (s_verbosity >= 3)
      pout() << "AMRLevelAdvectDiffuse::poissonSolveImplicitComposite1 " << m_level << endl;
      
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    
    int finest_level = hierarchy.size()-1;
    
    // Vector of pointers to LevelData of FABS
    Vector<LevelData<FArrayBox>* > phi(finest_level+1, NULL);
    Vector<LevelData<FArrayBox>* > rhs(finest_level+1, NULL);
    // srs is to store everything except the first term on the RHS of eq. (4.68) in Colella et al. (1999b)
    Vector<LevelData<FArrayBox>* > srs(finest_level+1, NULL);
    // dge is to store the negative divergence of deltaGamma_e of 4.67 and 4.68 in Colella et al. (1999b), which is used also to update density; note the stored divergence is multiplied by m_dt
    Vector<LevelData<FArrayBox>* > dge(finest_level+1, NULL);
    
    // Edge field correction
    Vector<LevelData<FluxBox>* > eec(finest_level+1, NULL);

    // loop over levels, allocate storage, set up for AMRSolve
    // if coarser level exists, define it as well for BCs.
    int startLev = Max(m_level-1, 0);
    
    for (int lev = finest_level; lev>= m_level; lev--)
      hierarchy[lev]->setPoissonCoeffABComposite1(m_dt);
    
    for (int lev = finest_level; lev>= startLev; lev--) {
      // rhs has no ghost cells, phi does
      rhs[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      phi[lev]  = new LevelData<FArrayBox>(grids[lev], 1, m_numGhost*IntVect::Unit);
      srs[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      dge[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      eec[lev]  = new LevelData<FluxBox>(grids[lev], 1, m_numGhost*IntVect::Unit);
      
      for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit) {
        (*phi[lev])[dit()].setVal(0.0);
        (*rhs[lev])[dit()].setVal(0.0);
        (*srs[lev])[dit()].setVal(0.0);
        (*dge[lev])[dit()].setVal(0.0);
        (*eec[lev])[dit()].setVal(0.0);
      }
    }
    
    for (int lev=m_level; lev>= m_level; lev--) {
      if (hierarchy[lev]->m_hasFiner && lev >= m_level) {
        AMRLevelAdvectDiffuse* cl = hierarchy[lev];
        AMRLevelAdvectDiffuse* fl = hierarchy[lev+1];
        
        // revisit: For cases with >2 levels, is dt here the one for lev or lbase?
        // Define coe = dt*neMidStep*mu
        LevelData<FluxBox> coe;
        Real scale = m_dt/cl->m_dt;
//        Real scale = 1.0;
        coe.define(*s_bCoef[lev]);
        (*s_bCoef[lev]).copyTo(coe);
        for (DataIterator dit=(*cl).m_grids.dataIterator(); dit.ok(); ++dit)
          coe[dit()] -= 1.0;
        
        cl->getDivDeltaFlux((*dge[lev]), &(cl->m_EEdgeOld), &(fl->m_EEdgeCoarseOld), &coe, -1.0/cl->m_dx * scale);
        cl->m_fluxRegister.reflux((*dge[lev]), -1.0/cl->m_dx * scale);
        // get average source from the next finer level
        fl->m_coarseAverage.averageToCoarse(*dge[lev], *dge[lev+1]);
        (*dge[lev]).copyTo(*srs[lev]);

        cl->getDivDeltaFlux((*srs[lev]), &(cl->m_EEdgeOld), &(fl->m_EEdgeCoarseOld), NULL, -1.0/cl->m_dx);
        cl->getDivDeltaFlux((*srs[lev]), &(cl->m_field.m_EEdge), &(fl->m_field.m_EEdge), NULL, 1.0/cl->m_dx);
        // get average source from the next finer level
        fl->m_coarseAverage.averageToCoarse(*srs[lev], *srs[lev+1]);
        // save the edge field from the predictor step for use by synchronization of next time step
        cl->m_field.m_EEdge.copyTo(cl->m_EEdgeOld);
        fl->m_field.m_EEdge.copyTo(fl->m_EEdgeCoarseOld);
      }
    }
    
//    if (m_level == 0 && s_verbosity >= 3) {
//      string filename("PoiImSrs");
//      filename.append(to_string(AMR::s_step+1));
//      filename.append(".hdf5");
//      WriteAMRHierarchyHDF5(filename, grids, srs, lev0Domain.domainBox(), refRat, hierarchy.size());
//    }
//    if (m_level == 0 && s_verbosity >= 3) {
//      string filename("dge");
//      filename.append(to_string(AMR::s_step+1));
//      filename.append(".hdf5");
//      WriteAMRHierarchyHDF5(filename, grids, dge, lev0Domain.domainBox(), refRat, hierarchy.size());
//    }
    
// assuming zero correction for the coarser level, which is used for BCs. That is no interpolation is done for m_level-1, if it exists.
    
    int numSmooth, numMG, maxIter, mgverb;
    Real tolerance, hang, normThresh;
    
    ParmParse ppPoisson("PoissonSolver");
    ppPoisson.get("num_smooth", numSmooth);
    ppPoisson.get("num_mg",     numMG);
    ppPoisson.get("hang_eps",   hang);
    ppPoisson.get("norm_thresh",normThresh);
    ppPoisson.get("tolerance",  tolerance);
    ppPoisson.get("max_iter",   maxIter);
    ppPoisson.get("verbosity",  mgverb);
    
    VCAMRPoissonOp2Factory* amrpopEPotVC = new VCAMRPoissonOp2Factory;
    // the BC function supplied below isn't used as homogeneous BCs are forced at the solve call.
    amrpopEPotVC->define(lev0Domain, grids, refRat, lev0Dx, m_EPotbcFunc, 0.0, s_aCoefComp, -1.0, s_bCoefComp);
    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > EPotImpOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(amrpopEPotVC);
    RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > EPotImpAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >
    (new AMRMultiGrid<LevelData<FArrayBox> >());
    EPotImpAMRMG->define(lev0Domain, *EPotImpOpFact, &s_EPotImpCompBotSolver, hierarchy.size());
    EPotImpAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
                                        maxIter, tolerance, hang, normThresh);
    EPotImpAMRMG->m_verbosity = mgverb;
    s_EPotImpCompBotSolver.m_verbosity = mgverb;
      
    int lbase = m_level;
    int lmax  = finest_level;;
    
    // for testing purpose, copy the field before synchronization for ouput
//    for (int lev = finest_level; lev >= m_level; lev--) {
//      hierarchy[lev]->m_field.copyTo(hierarchy[lev]->m_fieldOld);
//    }
    
    Real totVolCharge = 0;
    Real chargeTol = 1e-10;
    Real surfChargeRelTol = 1e-4;
    Real totSurfCharge = 0;
    int iter = 0;
    int maxNewtonIter = 5;
    Real relChange = 1;
    Real relChangeTol = 1e-3;
    
    {
      // use rhs as a temporary variable to store charge
      for (int lev = lbase; lev<= lmax; lev++) {
        for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
          (*rhs[lev])[dit()].setVal(0.0);
                
        hierarchy[lev]->m_ionNew.copyTo(Interval(0, 0), *rhs[lev], Interval(0, 0));
        for (DataIterator dit=rhs[lev]->dataIterator(); dit.ok(); ++dit) {
          (*rhs[lev])[dit()].minus(hierarchy[lev]->m_ionNew[dit()], 1, 0);
          (*rhs[lev])[dit()] -= hierarchy[lev]->m_UNew[dit()];
        }
      }
      // find the total volume charge on the finer levels
      totVolCharge = computeNorm(rhs, refRat, hierarchy[lbase+1]->m_dx, Interval(0, 0), 1, lbase+1);
      
      for (int lev = lbase; lev<= lmax; lev++) {
        for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
          (*rhs[lev])[dit()].setVal(0.0);
        if (hierarchy[lev]->m_hasFiner)
          hierarchy[lev]->getDivDeltaFlux(*rhs[lev], &(hierarchy[lev]->m_field.m_EEdge), &(hierarchy[lev+1]->m_field.m_EEdge), NULL, 1.0/hierarchy[lev]->m_dx);
      }
      totSurfCharge = computeNorm(rhs, refRat, hierarchy[lbase]->m_dx, Interval(0, 0), 1, lbase);
      
  //    if (m_level == 0 && s_verbosity >= 3) {
  //      string filename("PoiImSrs");
  ////      filename.append(to_string(AMR::s_step+1));
  //      filename.append(to_string(iter));
  //      filename.append(".hdf5");
  //      WriteAMRHierarchyHDF5(filename, grids, rhs, lev0Domain.domainBox(), refRat, hierarchy.size());
  //    }
    }
    if (s_verbosity >= 3) {
      pout() << "lbase = " << lbase << endl;
      pout() << "iter = " << iter << " " << "totVolCharge = " << totVolCharge << "  " << "totSurfCharge = " << totSurfCharge << endl;
    }
    // save the predictor solution
    for (int lev = m_level; lev <= finest_level; lev++)
      hierarchy[lev]->m_phi.copyTo(hierarchy[lev]->m_phiOld);
    
    if (m_EPotBCVarying)
      for (int lev = lbase; lev <= lmax; lev++) {
        // all these levels are at the same time; no need to use individual levels' times
        setTimeHelper(hierarchy[lev]->m_EPotbcFunc, m_time+m_dt);
        // pout() << "m_time+m_dt = " << m_time+m_dt << endl;
      }
    
//    {
//      std::ostringstream oss_phi;
//      oss_phi << "results/phiSync_lstep" << AMR::s_step
//              << "_lbase" << lbase
//              << "_iter" << iter << ".h5";
//      WriteAMRHierarchyHDF5(oss_phi.str(), grids, phi, lev0Domain.domainBox(), refRat, hierarchy.size());
//
//      std::ostringstream oss_rhs;
//      oss_rhs << "results/rhs_lstep" << AMR::s_step
//              << "_lbase" << lbase
//              << "_iter" << iter << ".h5";
//      WriteAMRHierarchyHDF5(oss_rhs.str(), grids, rhs, lev0Domain.domainBox(), refRat, hierarchy.size());
//    }


    while ((totSurfCharge > surfChargeRelTol * (totVolCharge+chargeTol) || iter == 0) && iter <= maxNewtonIter && relChange > relChangeTol) {
      // note rhs contains additional terms that depend on the correction
      for (int lev = finest_level; lev>= m_level; lev--) {
        (*srs[lev]).copyTo((*rhs[lev]));
        if (hierarchy[lev]->m_hasFiner) {
          hierarchy[lev]->getDivDeltaFlux(*rhs[lev], eec[lev], eec[lev+1], s_bCoefComp[lev], 1.0/hierarchy[lev]->m_dx);
          // get average source from the next finer level
          hierarchy[lev+1]->m_coarseAverage.averageToCoarse(*rhs[lev], *rhs[lev+1]);
        }
      }

      EPotImpAMRMG->solve(phi, rhs, lmax, lbase, false, true);
      
      for (int lev = m_level; lev <= finest_level; lev++) {
        hierarchy[lev]->m_phiOld.copyTo(hierarchy[lev]->m_phi);
        for (DataIterator dit=phi[lev]->dataIterator(); dit.ok(); ++dit) {
          hierarchy[lev]->m_phi[dit()] += (*phi[lev])[dit()];
        }
      }

      Real eps = 1.0e-10;
      { // find fields
        // find the edge field correction; assuming the coarser level of current level has phiCorr = 0;
        for (int lev= m_level; lev <= finest_level; lev++) {
          if (hierarchy[lev]->m_hasCoarser)
            hierarchy[lev]->fillGhost(hierarchy[lev]->m_pwl, *phi[lev], *phi[lev-1], 0, *phi[lev-1], 0, eps, 0);
          else
            (*phi[lev]).exchange();
          
          if (hierarchy[lev]->m_hasCoarser)
            Gradient::levelGradientMAC(*eec[lev], *phi[lev], phi[lev-1], hierarchy[lev]->m_dx, refRat[lev-1], hierarchy[lev]->m_problem_domain);
          else
            Gradient::levelGradientMAC(*eec[lev], *phi[lev], NULL, hierarchy[lev]->m_dx, 0, hierarchy[lev]->m_problem_domain);
          for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit) {
            (*eec[lev])[dit()].negate();
          }
        }

        // multigrid smoothes down first and then up; fine level solution is more accurate! This is really necessary; otherwise, there would be discontinuity at the coarse-fine interface!
        for (int lev = finest_level; lev >= m_level; lev--)
          if (hierarchy[lev]->m_hasCoarser)
            hierarchy[lev]->m_coarseAverage.averageToCoarse(hierarchy[lev-1]->m_phi, hierarchy[lev]->m_phi);
        
        Real time_eps = 1.0e-10;
        if (m_coarser_level_ptr)
          if (abs(m_coarser_level_ptr->time() - m_time) > time_eps)
            computeEField(true);
          else
            computeEField(false);
        else
          computeEField(false);
        
        for (int lev = m_level+1; lev <= finest_level; lev++)
          hierarchy[lev]->computeEField(false);
      }
        
      for (int lev = lbase; lev<= lmax; lev++) {
        for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
          (*rhs[lev])[dit()].setVal(0.0);
        if (hierarchy[lev]->m_hasFiner)
          hierarchy[lev]->getDivDeltaFlux(*rhs[lev], &(hierarchy[lev]->m_field.m_EEdge), &(hierarchy[lev+1]->m_field.m_EEdge), NULL, 1.0/hierarchy[lev]->m_dx);
      }
      relChange = totSurfCharge;
      totSurfCharge = computeNorm(rhs, refRat, hierarchy[lbase]->m_dx, Interval(0, 0), 1, lbase);
      relChange = fabs((totSurfCharge-relChange)/relChange);
      
      iter++;
      
//      {
//        std::ostringstream oss_phi;
//        oss_phi << "results/phiSync_step" << AMR::s_step
//                << "_lbase" << lbase
//                << "_iter" << iter << ".h5";
//        WriteAMRHierarchyHDF5(oss_phi.str(), grids, phi, lev0Domain.domainBox(), refRat, hierarchy.size());
//
//        std::ostringstream oss_rhs;
//        oss_rhs << "results/rhs_step" << AMR::s_step
//                << "_lbase" << lbase
//                << "_iter" << iter << ".h5";
//        WriteAMRHierarchyHDF5(oss_rhs.str(), grids, rhs, lev0Domain.domainBox(), refRat, hierarchy.size());
//      }

      
//      if (m_level == 0 && s_verbosity >= 3) {
//        string filename("PoiImSrs");
//  //      filename.append(to_string(AMR::s_step+1));
//        filename.append(to_string(iter));
//        filename.append(".hdf5");
//        WriteAMRHierarchyHDF5(filename, grids, rhs, lev0Domain.domainBox(), refRat, hierarchy.size());
//      }
            
      if (s_verbosity >= 3)
        pout() << "iter = " << iter << " " << "totVolCharge = " << totVolCharge << "  " << "totSurfCharge = " << totSurfCharge << "  " << "relChange = " << relChange << endl;
    }
    
//    for (int lev = finest_level; lev >= m_level; lev--) {
//      EdgeToCell(*eec[lev], hierarchy[lev]->m_fieldOld.m_E);
//    }
    
    // update density
    for (int lev = m_level; lev <= finest_level; lev++) {
      AMRLevelAdvectDiffuse* cl = hierarchy[lev];
      LevelData<FluxBox> fluxCorr;
      LevelData<FArrayBox> dU;
      fluxCorr.define(cl->m_grids, 1, IntVect::Zero);
      dU.define(cl->m_grids, 1, IntVect::Zero);
      (*s_bCoefComp[lev]).copyTo(fluxCorr);
      for (DataIterator dit=cl->m_grids.dataIterator(); dit.ok(); ++dit) {
        Box curBox = cl->m_grids.get(dit());
        FluxBox& curFlux = fluxCorr[dit()];
        FluxBox& curEEdge = (*eec[lev])[dit()];
        
        // bCoefComp contains additional 1.0, and dt at the coarsest level is already included.
        curFlux -= 1.0;
        curFlux *= curEEdge;
        curFlux.negate();
        FArrayBox& curdU = dU[dit()];
        curdU.setVal(0.0);
        for (int idir = 0; idir < SpaceDim; idir++) {
          // Compute flux difference fHi-fLo
          FArrayBox diff(curBox, dU.nComp());
          diff.setVal(0.0);

          FORT_FLUXDIFFF(CHF_FRA(diff),
                         CHF_CONST_FRA(curFlux[idir]),
                         CHF_CONST_INT(idir),
                         CHF_BOX(curBox));

          // Add flux difference to dU
          curdU += diff;
        }
        curdU *= -1.0/cl->m_dx;
        cl->m_UNew[dit()] += curdU;
        cl->m_UNew[dit()] += (*dge[lev])[dit()];
        // TODO: double check the following update to electron flux for current calculation
        curFlux *= (1/m_dt);
        cl->m_J.m_EEdge[dit()] += curFlux;
      }
      
      if (s_verbosity >= 3) {
        printDiagnosticInfo (cl->m_level, cl->m_dx, cl->m_grids, cl->m_UNew, "UNew", "postTimeStep");
        printDiagnosticInfo (cl->m_level, cl->m_dx, cl->m_grids, dU, "dU", "postTimeStep");
        printDiagnosticInfo (cl->m_level, cl->m_dx, cl->m_grids, fluxCorr, "flux", "postTimeStep");
        printDiagnosticInfo (cl->m_level, cl->m_dx, cl->m_grids, *eec[lev], "EEdge", "postTimeStep");
      }
    }
    
    //remember that startLev can be different from m_level
    for (int lev = startLev; lev<= finest_level; lev++) {
      delete rhs[lev];
      delete phi[lev];
      delete srs[lev];
      delete dge[lev];
      delete eec[lev];
      rhs[lev] = NULL;
      phi[lev] = NULL;
      srs[lev] = NULL;
      dge[lev] = NULL;
      eec[lev] = NULL;
    }
  }
}


void
AMRLevelAdvectDiffuse::
setPoissonCoeffABComposite1(Real commonDt) {
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::setPoissonCoeffABComposite1 " << m_level << endl;
  // on the edges contained in a c-f interface, dt of the coarser level should be used for calculating b coefficient
  (*s_aCoef[m_level]).copyTo(*s_aCoefComp[m_level]);
  (*s_bCoef[m_level]).copyTo(*s_bCoefComp[m_level]);
  
  Real eps = 1e-10;
  if (abs(m_dt-commonDt) > eps)
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
      (*s_bCoefComp[m_level])[dit()] -= 1.0;
      (*s_bCoefComp[m_level])[dit()] *= 1/m_dt;
      (*s_bCoefComp[m_level])[dit()] *= commonDt;
      (*s_bCoefComp[m_level])[dit()] += 1;
    }
  
  // on each edge in a c-f interface, the bCoeff on the coarser level the average of the finer level
  if (m_hasFiner) {
    AMRLevelAdvectDiffuse* amrGodFinerPtr = getFinerLevel();
//    (*s_bCoefComp[m_level+1]).exchange();
    amrGodFinerPtr->m_coarseAverageFace.averageToCoarse(*s_bCoefComp[m_level], *s_bCoefComp[m_level+1]);

  }
  
  if (s_verbosity >= 3) {
     printDiagnosticInfo (m_level, m_dx, m_grids, *s_bCoefComp[m_level], "bCoeffComp", "setPoissonCoeffABComp1");
  }
}
