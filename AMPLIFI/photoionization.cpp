//
//  photoionization.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#include <iostream>
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"
#include "UsingNamespace.H"

#include "additionalBCFunc.hpp"
#include "photoionization.hpp"

int photoionization::s_verbosity;
vector<double> SP3A{6.7e-3*1e2, 0.0346*1e2, 0.3059*1e2};
vector<double> SP3Lambda{0.0447*1e2, 0.1121*1e2, 0.5994*1e2};

void
PISetupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver) {
  
  CH_TIME("PISetupSolver");
  
  ParmParse ppSolver("PISolver");
  
  // multigrid solver parameters
  int numSmooth, numMG, maxIter, verbosity;
  Real eps, hang, normThresh = 1.0e-30;
  ppSolver.get("num_smooth",    numSmooth);
  ppSolver.get("tolerance",     eps);
  ppSolver.get("num_mg",        numMG);
  ppSolver.get("norm_thresh",   normThresh);
  ppSolver.get("hang_eps",      hang);
  ppSolver.get("max_iter",      maxIter);
  ppSolver.get("verbosity",     verbosity);
  
  a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                   numMG, maxIter, eps, hang, normThresh);
  a_amrSolver->m_verbosity = verbosity;
  
  // optional parameters
  ppSolver.query("num_pre", a_amrSolver->m_pre);
  ppSolver.query("num_post", a_amrSolver->m_post);
  ppSolver.query("num_bottom", a_amrSolver->m_bottom);
}

vector<int>    PIBCLo = vector<int>(SpaceDim, 0);
vector<int>    PIBCHi = vector<int>(SpaceDim, 0);
double         PIDiriNeumBCValues;

vector<double> PIRobinBCValues(3, 0);

void PIParseValue(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values) {
  a_values[0]=PIDiriNeumBCValues;
}

void PIParseValueRobin(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values) {
  a_values[0]=PIRobinBCValues[*side]; //store alpha
  a_values[1]=PIRobinBCValues[2];     //store gamma
}

void PIParseBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous) {
  
  if (!a_domain.domainBox().contains(a_state.box())) {
    
    Box valid = a_valid;
    int a_order = 1;
    
    for (int i=0; i<CH_SPACEDIM; ++i)
      
      // don't do anything if periodic
      if (!a_domain.isPeriodic(i)) {
      
        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        
        if (!a_domain.domainBox().contains(ghostBoxLo))
          
          switch (PIBCLo[i]) {
            case 0:
              DiriBC(a_state, valid, a_dx, a_homogeneous, PIParseValue, i, Side::Lo);
              break;
              
            case 1:
              NeumBC(a_state, valid, a_dx, a_homogeneous, PIParseValue, i, Side::Lo);
              break;
              
            case 2:
              ExtrapolateBC(a_state, valid, a_dx, i, Side::Lo, a_order);
              break;
            
            case 3:
              RobinBC(a_state, valid, a_dx, a_homogeneous, PIParseValueRobin, i, Side::Lo);
              break;
              
            default:
              MayDay::Error("bogus bc flag lo");
              break;
          }
        
        if (!a_domain.domainBox().contains(ghostBoxHi))
          
          switch (PIBCHi[i]) {
            case 0:
              DiriBC(a_state, valid, a_dx, a_homogeneous, PIParseValue, i, Side::Hi);
              break;
              
            case 1:
              NeumBC(a_state, valid, a_dx, a_homogeneous, PIParseValue, i, Side::Hi);
              break;
              
            case 2:
              ExtrapolateBC(a_state, valid, a_dx, i, Side::Hi, a_order);
              break;
              
            case 3:
              RobinBC(a_state, valid, a_dx, a_homogeneous, PIParseValueRobin, i, Side::Hi);
            break;
              
            default:
              MayDay::Error("bogus bc flag hi");
              break;
          }
        
      }
    // end if is not periodic in ith direction
  }
}


void photoionization::define(DisjointBoxLayout a_grids, int a_ncomps, IntVect ivGhost) {
  ParmParse pp;
  pp.get("verbosity",  s_verbosity);
  ncomps = a_ncomps;
  Psi.define(a_grids, ncomps, ivGhost);
  rate.define(a_grids, 1, ivGhost);
}

void photoionization::define(const photoionization& a_PI) {
  DisjointBoxLayout grids = a_PI.Psi.disjointBoxLayout();
  IntVect ivGhost = a_PI.Psi.ghostVect();
  define(grids, a_PI.ncomps, ivGhost);
}

void photoionization::copyTo(photoionization& des) {
  Psi.copyTo(Psi.interval(), des.Psi, des.Psi.interval());
  rate.copyTo(rate.interval(), des.rate, des.rate.interval());
  des.setCoefficients(A, lambda);
}

void photoionization::calcRate() {
  const DisjointBoxLayout& levelGrids = rate.getBoxes();
  DataIterator levelDit = levelGrids.dataIterator();
  for (levelDit.begin(); levelDit.ok(); ++levelDit) {
    FArrayBox& r = rate[levelDit];
    FArrayBox& p = Psi[levelDit];
    r.setVal(0.0);
    for (int comp = 0; comp < ncomps; comp++)
      r.plus(p, A[comp]*(constants::c0/lBar*tBar), comp, 0, 1);
  }
}

void photoionization::setCoefficients(vector<double> a_A, vector<double> a_lambda){
  A = a_A;
  lambda = a_lambda;
}

void PISetTestRHS(Vector<LevelData<FArrayBox>* > a_rhs, Vector<int>& a_refRat, ProblemDomain a_lev0Dom, Real a_lev0Dx, int startLev, int a_finestLevel, PITestProbType s_probtype) {
  
  CH_TIME("setRHS");
  
  Real levDx = a_lev0Dx;
  
  for (int lev=0; lev<startLev; lev++) {
    levDx /= a_refRat[lev];
  }
  
  for (int lev=startLev; lev<=a_finestLevel; lev++) {
    LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
    const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();
    // rhs is cell-centered...
    RealVect ccOffset = 0.5*levDx*RealVect::Unit;
    
    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit) {
      FArrayBox& thisRhs = levelRhs[levelDit];
      
      if (s_probtype == PITestProbType::zeroRHS)
        thisRhs.setVal(0.0);
      else if (s_probtype == PITestProbType::unityRHS)
        thisRhs.setVal(1.0);
      else if (s_probtype == PITestProbType::gaussian) {
        int numgaussian = 1;
        Vector<RealVect> center(numgaussian,RealVect::Zero);
        Vector<Real> scale(numgaussian, 1.0);
        Vector<Real> strength(numgaussian, 1.0);
        
        strength[0] = 3.5e22*1e6*pow(lBar, 3)*tBar;
        scale[0] = 0.01e-2/lBar;
        center[0] = 0.5*a_lev0Dx*RealVect(a_lev0Dom.size());
        
        thisRhs.setVal(0.0);
        
        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit) {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= levDx;
          loc += ccOffset;
          
          for (int n=0; n<numgaussian; n++) {
            RealVect dist = loc - center[n];
            Real radSqr = D_TERM(dist[0]*dist[0],
                                 +dist[1]*dist[1],
                                 +dist[2]*dist[2]);
            
            Real val = strength[n]*exp(-radSqr/(scale[n]*scale[n]));
            thisRhs(iv,0) += val;
          }
        }
      }
      else
        MayDay::Error("undefined problem type");
    } // end loop over grids on this level
    levDx /= a_refRat[lev];
  } // end loop over levels
}

