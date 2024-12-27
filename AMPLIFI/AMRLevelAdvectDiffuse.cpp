// code adapted from the AdvectDiffuse example in the Chombo package

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

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

#include "additionalBCFunc.hpp"
#include "photoionization.hpp"
#include "physicalConstants.h"
#include "Gradient.H"
#include "dataFileIFReduced.hpp"

RefCountedPtr<LevelTGA>                                   AMRLevelAdvectDiffuse::s_diffuseLevTGA = RefCountedPtr<LevelTGA>();
RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >       AMRLevelAdvectDiffuse::s_diffuseAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_diffuseOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
BiCGStabSolver<LevelData< FArrayBox> >                    AMRLevelAdvectDiffuse::s_botSolver;

RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >       AMRLevelAdvectDiffuse::s_EPotAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_EPotOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
BiCGStabSolver<LevelData< FArrayBox> >                    AMRLevelAdvectDiffuse::s_EPotBotSolver;

//RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >       AMRLevelAdvectDiffuse::s_EPotImpAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
//RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_EPotImpOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
BiCGStabSolver<LevelData< FArrayBox> >                    AMRLevelAdvectDiffuse::s_EPotImpBotSolver;
BiCGStabSolver<LevelData< FArrayBox> >                    AMRLevelAdvectDiffuse::s_EPotImpCompBotSolver;

RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >       AMRLevelAdvectDiffuse::s_PIAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_PIOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
BiCGStabSolver<LevelData< FArrayBox> >                    AMRLevelAdvectDiffuse::s_PIBotSolver;

Vector<RefCountedPtr<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_aCoef(10);
Vector<RefCountedPtr<LevelData<FluxBox> > >  AMRLevelAdvectDiffuse::s_bCoef(10);
Vector<RefCountedPtr<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_aCoefComp(10);
Vector<RefCountedPtr<LevelData<FluxBox> > >  AMRLevelAdvectDiffuse::s_bCoefComp(10);

bool AMRLevelAdvectDiffuse::s_testing;
unsigned long long AMRLevelAdvectDiffuse::totNumAdvCells = 0;
/*******/
AMRLevelAdvectDiffuse::
~AMRLevelAdvectDiffuse()
{
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::~AMRLevelAdvectDiffuse " << m_level << endl;
  s_diffuseLevTGA  = RefCountedPtr<LevelTGA>();
  s_diffuseAMRMG   = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
  s_diffuseOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
  s_EPotAMRMG      = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > ();
  s_EPotOpFact     = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > ();
  if(m_doImplicitPoisson) {
//    s_EPotImpAMRMG      = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > ();
//    s_EPotImpOpFact     = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > ();
    s_aCoef.pop_back();
    s_bCoef.pop_back();
    s_aCoefComp.pop_back();
    s_bCoefComp.pop_back();
  }
  s_PIAMRMG      = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > ();
  s_PIOpFact     = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > ();
  // assuming the last level should be destructed
}

/*******/
void
AMRLevelAdvectDiffuse::
getHierarchyAndGrids(Vector<AMRLevelAdvectDiffuse*>&        a_hierarchy,
                     Vector<DisjointBoxLayout>&             a_grids,
                     Vector<int>&                           a_refRat,
                     ProblemDomain&                         a_lev0Dom,
                     Real&                                  a_lev0Dx)
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int nlevels = hierarchy.size();
  
  a_hierarchy.resize(nlevels);
  a_refRat.resize(nlevels);
  a_grids.resize(nlevels);
  
  AMRLevelAdvectDiffuse* coarsestLevel = (AMRLevelAdvectDiffuse*)(hierarchy[0]);
  a_lev0Dx       = coarsestLevel->m_dx;
  a_lev0Dom      = coarsestLevel->m_problem_domain;
  
  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    AMRLevelAdvectDiffuse* adLevel = (AMRLevelAdvectDiffuse*)(hierarchy[ilev]);
    
    a_hierarchy[ilev] = adLevel;
    a_grids [ilev] = adLevel->m_grids;
    a_refRat[ilev] = adLevel->m_ref_ratio;
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
defineSolvers()
{
  
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::defineSolvers " << m_level << endl;

  int numSmooth, numMG, maxIter, mgverb;
  Real tolerance, hang, normThresh;
  
  ParmParse pp("diffusionSolver");
  pp.get("num_smooth", numSmooth);
  pp.get("num_mg",     numMG);
  pp.get("hang_eps",   hang);
  pp.get("norm_thresh",normThresh);
  pp.get("tolerance",  tolerance);
  pp.get("max_iter",   maxIter);
  pp.get("verbosity",  mgverb);
  
  Vector<AMRLevelAdvectDiffuse*>  hierarchy;
  Vector<DisjointBoxLayout>       grids;
  Vector<int>                     refRat;
  ProblemDomain                   lev0Dom;
  Real                            lev0Dx;
  getHierarchyAndGrids(hierarchy, grids, refRat, lev0Dom, lev0Dx);
  
  if (m_hasDiffusion) {
    s_botSolver.m_verbosity = mgverb-3;
    
    AMRPoissonOpFactory* amrpop =  new AMRPoissonOpFactory();
    //alpha = 1.0; See Chombo 4.2
    amrpop->define(lev0Dom, grids, refRat, lev0Dx, m_bcFunc, 1.0, m_nu);
    s_diffuseOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(amrpop);
    
    s_diffuseAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >
    (new AMRMultiGrid<LevelData<FArrayBox> >());
    
    s_diffuseAMRMG->define(lev0Dom, *s_diffuseOpFact, &s_botSolver, hierarchy.size());
    s_diffuseAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
                                        maxIter, tolerance, hang, normThresh);
    s_diffuseAMRMG->m_verbosity = mgverb;
    
    s_diffuseLevTGA = RefCountedPtr<LevelTGA>
    (new LevelTGA(grids, refRat, lev0Dom, s_diffuseOpFact, s_diffuseAMRMG));
  }
  
  ParmParse ppPoisson("PoissonSolver");
  ppPoisson.get("num_smooth", numSmooth);
  ppPoisson.get("num_mg",     numMG);
  ppPoisson.get("hang_eps",   hang);
  ppPoisson.get("norm_thresh",normThresh);
  ppPoisson.get("tolerance",  tolerance);
  ppPoisson.get("max_iter",   maxIter);
  ppPoisson.get("verbosity",  mgverb);
  
  s_EPotBotSolver.m_verbosity = mgverb;
  
  AMRPoissonOpFactory* amrpopEPot =  new AMRPoissonOpFactory();
  amrpopEPot->define(lev0Dom, grids, refRat, lev0Dx, m_EPotbcFunc, 0.0, 1.0);
  s_EPotOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(amrpopEPot);
  s_EPotAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >
  (new AMRMultiGrid<LevelData<FArrayBox> >());
  s_EPotAMRMG->define(lev0Dom, *s_EPotOpFact, &s_EPotBotSolver, hierarchy.size());
  s_EPotAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
                                   maxIter, tolerance, hang, normThresh);
  s_EPotAMRMG->m_verbosity = mgverb;
  
  s_PIBotSolver.m_verbosity = mgverb;
  AMRPoissonOpFactory* amrpopPI =  new AMRPoissonOpFactory();
  amrpopPI->define(lev0Dom, grids, refRat, lev0Dx, m_phtznbcFunc, 1.0, 1.0);
  s_PIOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(amrpopPI);
  s_PIAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >
  (new AMRMultiGrid<LevelData<FArrayBox> >());
  s_PIAMRMG->define(lev0Dom, *s_PIOpFact, &s_PIBotSolver, hierarchy.size());
  PISetupSolver(s_PIAMRMG);
  
  ParmParse pph("PISolver");
  pph.get("runSolve", m_phtzn.runSolve);
}

/********/
void
AMRLevelAdvectDiffuse::
define(const AdvectPhysics&        a_gphys,
       gas                         a_gas,
       BCHolder                    a_bcFunc,
       BCHolder                    a_EPotbcFunc,
       BCHolder                    a_PIbcFunc,
       const Real&                 a_cfl,
       const Real&                 a_domainLength,
       const Real&                 a_refineThresh,
       const int&                  a_tagBufferSize,
       const Real&                 a_initialDtMultiplier,
       const bool&                 a_useLimiting,
       const Real&                 a_nu)
{
  int a_verbosity;
  ParmParse pp;
  pp.get("verbosity",  a_verbosity);
  bool imReflux;
  pp.get("implicitReflux", imReflux);
  pp.get("sourceNumericalScheme", m_sourceNumericalScheme);
  verbosity(a_verbosity);
  pp.get("testing", s_testing);
  
  m_isDefined = true;
  m_cfl = a_cfl;
  m_domainLength = a_domainLength;
  m_refineThresh = a_refineThresh;
  m_tagBufferSize = a_tagBufferSize;
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initial_dt_multiplier = a_initialDtMultiplier;
  m_useLimiting = a_useLimiting;
  m_nu = a_nu;
  m_doImplicitReflux = imReflux;
  m_hasDiffusion     = (m_nu > 0);
   
  m_gas = a_gas;
  m_bcFunc  = a_bcFunc;
  ParmParse ppPoisson("PoissonSolver");
  if (ppPoisson.contains("varyingField"))
    ppPoisson.get("varyingField", m_varyingField);
  else
    m_varyingField = true;
  ppPoisson.get("implicit", m_doImplicitPoisson);
  m_EPotbcFunc = a_EPotbcFunc;
  m_phtznbcFunc = a_PIbcFunc;
  m_phtzn.setCoefficients(SP3A, SP3Lambda);
  GodunovPhysics* gphysPtr = a_gphys.new_godunovPhysics();
  m_advPhys = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtr);
}

/********/
void AMRLevelAdvectDiffuse::define(AMRLevel*            a_coarserLevelPtr,
                                   const ProblemDomain& a_problemDomain,
                                   int                  a_level,
                                   int                  a_refRatio)
{
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::define " << m_level << endl;
  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);
  
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelAdvectDiffuse* amrGodPtr = dynamic_cast<AMRLevelAdvectDiffuse*>(a_coarserLevelPtr);
    
    if (amrGodPtr != NULL)
    {
      define(*amrGodPtr->m_advPhys,
             amrGodPtr->m_gas,
             amrGodPtr->m_bcFunc,
             amrGodPtr->m_EPotbcFunc,
             amrGodPtr->m_phtznbcFunc,
             amrGodPtr->m_cfl,
             amrGodPtr->m_domainLength,
             amrGodPtr->m_refineThresh,
             amrGodPtr->m_tagBufferSize,
             amrGodPtr->m_initialDtMultiplier,
             amrGodPtr->m_useLimiting,
             amrGodPtr->m_nu);
    }
    else
    {
      MayDay::Error("AMRLevelAdvectDiffuse::define: a_coarserLevelPtr is not castable to AMRLevelAdvectDiffuse*");
    }
  }
  
  // Compute the grid spacing
  m_dx = m_domainLength/a_problemDomain.domainBox().longside();
  
  m_numGhost = 4;
  m_stateNames  = Vector<string>(1, string("ele"));
  m_advPhys->define(m_problem_domain, m_dx);
  PhysIBC* physIBCPtr = m_advPhys->getPhysIBC();
  physIBCPtr->define(m_problem_domain, m_dx);
}

/********/
void
AMRLevelAdvectDiffuse::
getCoarseDataPointers(LevelData<FArrayBox>** a_coarserDataOldPtr,
                      LevelData<FArrayBox>** a_coarserDataNewPtr,
                      LevelFluxRegister**    a_coarserFRPtr,
                      LevelFluxRegister**    a_finerFRPtr,
                      Real& a_tCoarserOld,
                      Real& a_tCoarserNew)
{
  *a_coarserDataOldPtr = NULL;
  *a_coarserDataNewPtr = NULL;
  *a_coarserFRPtr = NULL;
  *a_finerFRPtr   = NULL;
  
  a_tCoarserOld = 0.0;
  a_tCoarserNew = 0.0;
  
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelAdvectDiffuse* coarserPtr = getCoarserLevel();
    
    // Recall that my flux register goes between my level and the next
    // finer level
    *a_coarserFRPtr = &coarserPtr->m_fluxRegister;
    
    *a_coarserDataOldPtr = &coarserPtr->m_UOld;
    *a_coarserDataNewPtr = &coarserPtr->m_UNew;
    
    a_tCoarserNew = coarserPtr->m_time;
    a_tCoarserOld = a_tCoarserNew - coarserPtr->m_dt;
  }
  
  // A finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level
    *a_finerFRPtr = &m_fluxRegister;
  }
}

void
AMRLevelAdvectDiffuse::
fillGhost(PiecewiseLinearFillPatch&   a_pwl,
          LevelData<FArrayBox>&       a_U,
          const LevelData<FArrayBox>& a_UCoarseOld,
          const Real&                 a_TCoarseOld,
          const LevelData<FArrayBox>& a_UCoarseNew,
          const Real&                 a_TCoarseNew,
          const Real&                 a_dt,
          const Real&                 a_time)
{
  a_U.exchange();
  
  // Fill a_U's ghost cells using fillInterp
  if (m_hasCoarser) {
    Real time_eps = 1.0e-20;
    Real alpha;
    
    if (a_TCoarseNew - a_TCoarseOld > time_eps)
      // Fraction "a_time" falls between the old and the new coarse times
      alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);
    else
      alpha = 1.0;
    
    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
    Real eps = 0.04 * a_dt / m_ref_ratio;
    
    if (Abs(alpha) < eps)
      alpha = 0.0;
    
    if (Abs(1.0-alpha) < eps)
      alpha = 1.0;
    
    // Current time before old coarse time
    if (alpha < 0.0)
      MayDay::Error( "LevelAdvect::step: alpha < 0.0");
    
    // Current time after new coarse time
    if (alpha > 1.0)
      MayDay::Error( "LevelAdvect::step: alpha > 1.0");
    
    // Interpolate ghost cells from next coarser level using both space
    // and time interpolation
    a_pwl.fillInterp(a_U, a_UCoarseOld, a_UCoarseNew, alpha, 0, 0, a_U.nComp());
  }
}

void
AMRLevelAdvectDiffuse::
fillGhostFace(LevelData<FluxBox>&       a_U,
              const LevelData<FluxBox>& a_UCoarseOld,
              const Real&               a_TCoarseOld,
              const LevelData<FluxBox>& a_UCoarseNew,
              const Real&               a_TCoarseNew,
              const Real&               a_dt,
              const Real&               a_time)
{
  a_U.exchange();
  
  // Fill a_U's ghost cells using fillInterp
  if (m_hasCoarser) {
    Real time_eps = 1.0e-20;
    Real alpha;
    
    if (a_TCoarseNew - a_TCoarseOld > time_eps)
      // Fraction "a_time" falls between the old and the new coarse times
      alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);
    else
      alpha = 1.0;
    
    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
    Real eps = 0.04 * a_dt / m_ref_ratio;
    
    if (Abs(alpha) < eps)
      alpha = 0.0;
    
    if (Abs(1.0-alpha) < eps)
      alpha = 1.0;
    
    // Current time before old coarse time
    if (alpha < 0.0)
      MayDay::Error( "LevelAdvect::step: alpha < 0.0");
    
    // Current time after new coarse time
    if (alpha > 1.0)
      MayDay::Error( "LevelAdvect::step: alpha > 1.0");
    
    // Interpolate ghost cells from next coarser level using both space
    // and time interpolation
    m_pwlf.fillInterp(a_U, a_UCoarseOld, a_UCoarseNew, alpha, 0, 0, a_U.nComp());
  }
}

void
AMRLevelAdvectDiffuse::
fillGhostQuad(QuadCFInterp&               a_qcfi,
              LevelData<FArrayBox>&       a_U,
              const LevelData<FArrayBox>& a_UCoarseOld,
              const Real&                 a_TCoarseOld,
              const LevelData<FArrayBox>& a_UCoarseNew,
              const Real&                 a_TCoarseNew,
              const Real&                 a_time)
{
  a_U.exchange();
  
  // Fill a_U's ghost cells using fillInterp
  if (m_hasCoarser) {
    LevelData<FArrayBox> UCoarseInterp;
    UCoarseInterp.define(a_UCoarseNew);
    interpolateInTime(UCoarseInterp, a_UCoarseOld, a_UCoarseNew, a_time, a_TCoarseOld, a_TCoarseNew);
    // Interpolate ghost cells from next coarser level using both space
    // and time interpolation
    a_qcfi.coarseFineInterp(a_U, UCoarseInterp);
  }
}

/*******/
Real
AMRLevelAdvectDiffuse::
advance()
{
  if (s_verbosity >= 2)
    pout() << "AMRLevelAdvectDiffuse::advance " << m_level << endl;
  
  // Copy the new to the old
  m_UNew.copyTo(m_UNew.interval(), m_UOld, m_UOld.interval());
  m_ionNew.copyTo(m_ionNew.interval(), m_ionOld, m_ionOld.interval());
  m_phi.copyTo(m_phi.interval(), m_phiOld, m_phiOld.interval());
  m_field.copyTo(m_fieldOld);
  m_mu.copyTo(m_muOld);
  m_advVel.copyTo(m_advVelOld);
  
  LevelData<FArrayBox> diffusiveSrc(m_grids, 1, IntVect::Unit);
  makeDiffusiveSource(diffusiveSrc, m_UOld);
  
  printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "AMRLevelAdvectDiffuse::advance");
  printDiagnosticInfo (m_level, m_dx, m_grids, m_ionNew, "ion", "AMRLevelAdvectDiffuse::advance");
  
  Real newDt = diffusiveAdvance(diffusiveSrc);
  
  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_cfl * newDt;
  m_dtNew = returnDt;
  
  return returnDt;
}

need to modify ion densities too
/*********/
void
AMRLevelAdvectDiffuse::
updateWithReactionContribution (LevelData<FArrayBox>& U, LevelData<FArrayBox>& UIon, const LevelData<FArrayBox>& Emag, const LevelData<FArrayBox>& mu, const double dt) {
  
  for (DataIterator dit=U.dataIterator(); dit.ok(); ++dit) {
    FArrayBox dU, rateI, rateA;
    dU.define(U[dit()].box(), U[dit()].nComp());
    rateI.define(U[dit()].box(), U[dit()].nComp());
    rateA.define(U[dit()].box(), U[dit()].nComp());
    if (m_gas.m_uniformity) {
      getRate(rateI, Emag[dit()], mu[dit()], "ionization");
      getRate(rateA, Emag[dit()], mu[dit()], "attachment");
    } else {
      getRate(rateI, Emag[dit()], m_neut[dit()], mu[dit()], "ionization");
      getRate(rateA, Emag[dit()], m_neut[dit()], mu[dit()], "attachment");
    }
    
    dU.setVal(0.0);
    dU.axby(rateI, rateA, dt, -dt);
    dU *= U[dit()];
    U[dit()] += dU;
  }
  U.exchange();
}

/*********/
void
AMRLevelAdvectDiffuse::
updateWithReactionContribution (LevelData<FArrayBox>& U, const LevelData<FArrayBox>& Emag, const LevelData<FArrayBox>& mu, const LevelData<FArrayBox>& phtznRate, const Real dt) {
  
  for (DataIterator dit=U.dataIterator(); dit.ok(); ++dit) {
    FArrayBox dU, rateI, rateA;
    dU.define(U[dit()].box(), U[dit()].nComp());
    rateI.define(U[dit()].box(), U[dit()].nComp());
    rateA.define(U[dit()].box(), U[dit()].nComp());
    if (m_gas.m_uniformity) {
      getRate(rateI, Emag[dit()], mu[dit()], "ionization");
      getRate(rateA, Emag[dit()], mu[dit()], "attachment");
    } else {
      getRate(rateI, Emag[dit()], m_neut[dit()], mu[dit()], "ionization");
      getRate(rateA, Emag[dit()], m_neut[dit()], mu[dit()], "attachment");
    }
    
    dU.setVal(0.0);
    dU.axby(rateI, rateA, 1.0, -1.0);
    dU *= U[dit()];
    dU += phtznRate[dit()];
    dU *= dt;
    U[dit()] += dU;
  }
  U.exchange();
}

/*********/
void
AMRLevelAdvectDiffuse::
makeDiffusiveSource(LevelData<FArrayBox>& a_diffusiveSrc, const LevelData<FArrayBox>& U) {
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    a_diffusiveSrc[dit()].setVal(0);
  if (m_hasDiffusion) {
    RefCountedPtr<AMRPoissonOp> amrpop = RefCountedPtr<AMRPoissonOp>((AMRPoissonOp*)s_diffuseOpFact->AMRnewOp(m_problem_domain));
    LevelData<FArrayBox> zero;
    amrpop->create(zero, a_diffusiveSrc);
    amrpop->setToZero(zero);
    // see AMRPoissonOP.cpp: m_alpha = a_alpha * m_aCoef; m_beta  = a_beta  * m_bCoef;
    amrpop->setAlphaAndBeta(0., 1.0);
    if (m_level == 0)
      amrpop->residual(a_diffusiveSrc, U, zero, false);
    else {
      LevelFluxRegister* coarserFRPtr=NULL;
      LevelFluxRegister* finerFRPtr  =NULL;
      LevelData<FArrayBox>* coarserDataOldPtr = NULL;
      LevelData<FArrayBox>* coarserDataNewPtr = NULL;
      Real tCoarserOld, tCoarserNew;
      
      getCoarseDataPointers(&coarserDataOldPtr,
                            &coarserDataNewPtr,
                            &coarserFRPtr,
                            &finerFRPtr,
                            tCoarserOld, tCoarserNew);
      
      LevelData<FArrayBox> UCoarse(coarserDataOldPtr->disjointBoxLayout(),
                                   coarserDataOldPtr->nComp(),
                                   coarserDataOldPtr->ghostVect());
      
      Real interpTime = m_time + 0.5*m_dt;
      interpolateInTime(UCoarse, *coarserDataOldPtr, *coarserDataNewPtr,
                        interpTime, tCoarserOld, tCoarserNew);
      amrpop->AMRResidualNF(a_diffusiveSrc, U, UCoarse, zero, false);
    }
    
    amrpop->scale(a_diffusiveSrc, -1.0);
    
    /// Over the coarse-fine interface, the diffusive source is set to zero. At fine-fine interfaces, it is filled in by neighboring boxes.
    a_diffusiveSrc.exchange();
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
interpolateInTime(LevelData<FArrayBox>&          a_interp,
                  const LevelData<FArrayBox>&    a_old,
                  const LevelData<FArrayBox>&    a_new,
                  Real a_time, Real a_tOld, Real a_tNew)
{
  CH_assert(a_tNew >= a_tOld);
  //interp = alpha* unew + (1-alpha) uold
  Real alpha = 0;
  Real diff = a_tNew-a_tOld;
  if (diff > 0)
  {
    CH_assert(a_time >= a_tOld);
    alpha = (a_time - a_tOld)/diff;
  }
  for (DataIterator dit = a_interp.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox tempInte(a_interp[dit()].box(), a_interp.nComp());
    tempInte.copy(a_new[dit()]);
    a_interp[dit()].copy(a_old[dit()]);
    tempInte        *= alpha;     //temp has alpha*unew
    a_interp[dit()] *= 1.0-alpha; //interp has (1-alpha)uold
    a_interp[dit()] += tempInte;  //interp now has (1-alpha)uold + alpha*unew
  }
}

/*********/
Real
AMRLevelAdvectDiffuse::
diffusiveAdvance(LevelData<FArrayBox>& a_diffusiveSrc)
{
  LevelFluxRegister* coarserFRPtr=NULL;
  LevelFluxRegister* finerFRPtr  =NULL;
  LevelData<FArrayBox>* coarserDataOldPtr = NULL;
  LevelData<FArrayBox>* coarserDataNewPtr = NULL;
  Real tCoarserOld, tCoarserNew;
  
  getCoarseDataPointers(&coarserDataOldPtr,
                        &coarserDataNewPtr,
                        &coarserFRPtr,
                        &finerFRPtr,
                        tCoarserOld, tCoarserNew);
  
  LevelData<FArrayBox> srs, srsTmp, srsRea, srsReaTmp; // divide sources into general source and reaction term
  srs.define(m_UNew);
  srsTmp.define(m_UNew);
  
  // call poissonSolve in case that the source charge at the current level is modified after last Poisson solve
  if (m_varyingField) {
    poissonSolve();
    fillMobility(false);
    fillAdvectionVelocity(false);
  }
  
  
//  makeSource(m_UNew, 0.0);
  
  
  for (DataIterator dit=srs.dataIterator(); dit.ok(); ++dit) {
    FArrayBox rateI, rateA;
    rateI.define(srs[dit()].box(), srs[dit()].nComp());
    rateA.define(srs[dit()].box(), srs[dit()].nComp());
    if (m_gas.m_uniformity) {
      getRate(rateI, m_field.m_Emag[dit()], m_mu[dit()], "ionization");
      getRate(rateA, m_field.m_Emag[dit()], m_mu[dit()], "attachment");
    } else {
      getRate(rateI, m_field.m_Emag[dit()], m_neut[dit()], m_mu[dit()], "ionization");
      getRate(rateA, m_field.m_Emag[dit()], m_neut[dit()], m_mu[dit()], "attachment");
    }
    
    srs[dit()].setVal(0.0);
    srsTmp[dit()].setVal(0.0);
    srs[dit()].axby(rateI, rateA, 1.0, -1.0);
    srs[dit()] *= m_UNew[dit()];
    srs[dit()] += m_phtzn.rate[dit()];
  }
  srs.exchange();
  
  srs.copyTo(srsTmp);

  for (DataIterator dit=srs.dataIterator(); dit.ok(); ++dit)
    srsTmp[dit()] += a_diffusiveSrc[dit()];
  
  //  patchgodunov L273
  //  localSource *= 0.5 * a_dt;
  // *finerFRPtr is set to zero at the beginning but *coarserFRPtr is incremented from prior value
  Real newDt = m_levelGodunov.step(m_UNew,
                                   m_flux,
                                   m_advVel,
                                   srsTmp,
                                   *coarserDataOldPtr,
                                   tCoarserOld,
                                   *coarserDataNewPtr,
                                   tCoarserNew,
                                   m_time,
                                   m_dt);
  
  if (m_hasFiner)
    finerFRPtr->setToZero();
  
  Interval interv(0, 0);
  if (m_hasDiffusion) {
    //compute Du = unew-uold and put uold back into unew
    //so we can advance using leveltga
    m_UNew.copyTo(interv, m_dU,   interv);
    for (DataIterator dit=m_dU.dataIterator(); dit.ok(); ++dit) {
      m_dU[dit()] -= m_UOld[dit()];
      m_dU[dit()] /= m_dt;
      m_dU[dit()] += srs[dit()];
    }
    s_diffuseLevTGA->computeDiffusion(m_dUDiff, m_UOld, m_dU,
                                finerFRPtr, coarserFRPtr,
                                coarserDataOldPtr, coarserDataNewPtr,
                                m_time, tCoarserOld, tCoarserNew,
                                m_dt, m_level);
  }
  
  m_UOld.copyTo(interv, m_UNew, interv);
  
  if (m_doImplicitPoisson) {
    // m_flux is used in Poisson solve to set the variable coefficients
    if (m_varyingField)
      poissonSolveImplicit();
    for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
      Real eps = 1e-20;
      const Box& b = m_grids[dit()];
      FluxBox& curFlux = m_flux[dit()];
      FluxBox& advVel = m_advVel[dit()];
      advVel += eps;
      curFlux.divide(advVel, b, 0, 0);
    }
    fillAdvectionVelocity(true);
    for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
      const Box& b = m_grids[dit()];
      FluxBox& curFlux = m_flux[dit()];
      FluxBox& advVel = m_advVel[dit()];
      curFlux.mult(advVel, b, 0, 0);
    }
    // mobility is updated after the drift flux is computed so that the same drift flux is used in both Poisson's equation and the continuity equation
    fillMobility(true);
  }
    
  for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit)
    m_UNew[dit()].plus(m_dUDiff[dit()], m_dt);
  
  for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
    Interval UInterval(0, m_UNew.nComp()-1);
        
    Box curBox = m_grids.get(dit());
    FluxBox& curFlux = m_flux[dit()];
    m_dU[dit()].setVal(0.0);
    for (int idir = 0; idir < SpaceDim; idir++) {
      // Compute flux difference fHi-fLo
      FArrayBox diff(m_dU[dit()].box(), m_dU.nComp());
      diff.setVal(0.0);

      FORT_FLUXDIFFF(CHF_FRA(diff),
                     CHF_CONST_FRA(curFlux[idir]),
                     CHF_CONST_INT(idir),
                     CHF_BOX(curBox));

      // Add flux difference to dU
      m_dU[dit()] += diff;
    }
    m_dU[dit()] *= -m_dt/m_dx;
    m_UNew[dit()] += m_dU[dit()];
    
    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++) {
      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner)
        (*finerFRPtr).incrementCoarse(curFlux[idir],m_dt,dit(),
                                      UInterval,
                                      UInterval,idir);
      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser)
        (*coarserFRPtr).incrementFine(curFlux[idir],m_dt,dit(),
                                      UInterval,
                                      UInterval,idir);
    }
  }
  
  if (!m_doImplicitPoisson && m_varyingField)
    // Get the provisonal potential and mobility at the end of this step, which are used to fill ghost cells of the next finer level
    if (m_varyingField) {
      poissonSolve();
      fillMobility(true);
      fillAdvectionVelocity(true);
    }
  
  switch (m_sourceNumericalScheme) {
    
    // Dhali and Williams, 1987, but with the rate calculated using average field
    case 2: {
      double alpha = 0.5;
      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateI(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateA(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        
        ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
        muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
        if (m_gas.m_uniformity) {
          getRate(rateI, ETmp, muTemp, "ionization");
          getRate(rateA, ETmp, muTemp, "attachment");
        } else {
          getRate(rateI, ETmp, m_neut[dit()], muTemp, "ionization");
          getRate(rateA, ETmp, m_neut[dit()], muTemp, "attachment");
        }
        
        FArrayBox Uhs(m_UNew[dit()].box(), m_UNew.nComp()); // U at half step
        Uhs.axby(m_UOld[dit()], m_UNew[dit()], 1-alpha, alpha);
        srs[dit()].axby(rateI, rateA, 1.0, -1.0);
       
        srs[dit()].mult(m_UOld[dit()]);
        
        srs[dit()].plus(m_phtzn.rate[dit()]);
        Uhs.plus(srs[dit()], 0.5*m_dt);
        
        srs[dit()].copy(rateI);
        srs[dit()].mult(Uhs);
        srs[dit()].plus(m_phtzn.rate[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        
        srs[dit()].copy(rateA);
        srs[dit()].mult(Uhs);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], -m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 1);
      }
      
      break;
    }
     
    // differ with 2 only in Uhs, where source is calcuated using U mean
    case 3:
    {
      double alpha = 0.5;
      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateI(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateA(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        
        ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
        muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
        if (m_gas.m_uniformity) {
          getRate(rateI, ETmp, muTemp, "ionization");
          getRate(rateA, ETmp, muTemp, "attachment");
        } else {
          getRate(rateI, ETmp, m_neut[dit()], muTemp, "ionization");
          getRate(rateA, ETmp, m_neut[dit()], muTemp, "attachment");
        }
        
        FArrayBox Uhs(m_UNew[dit()].box(), m_UNew.nComp()); // U at half step
        Uhs.axby(m_UOld[dit()], m_UNew[dit()], 1-alpha, alpha);
        srs[dit()].axby(rateI, rateA, 1.0, -1.0);
        
        srs[dit()].mult(Uhs);
        
        srs[dit()].plus(m_phtzn.rate[dit()]);
        Uhs.plus(srs[dit()], 0.5*m_dt);
        
        srs[dit()].copy(rateI);
        srs[dit()].mult(Uhs);
        srs[dit()].plus(m_phtzn.rate[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        
        srs[dit()].copy(rateA);
        srs[dit()].mult(Uhs);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], -m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 1);
      }
      
      break;
    }
      
//    // LeVeque, p. 390, 2002
//    case 4: {
//      double rateMax = 0;
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        FArrayBox rateNew(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//        rateNew.setVal(0.0);
//        getRate(rateNew, m_field.m_Emag[dit()], m_mu[dit()], "ionization");
//        rateMax = max(rateMax, rateNew.max());
//      }
//
//      LevelData<FArrayBox> dUTmp;
//      dUTmp.define(m_UOld);
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        dUTmp[dit()].setVal(0.0);
//      }
//
//      int numOfIters = m_dt*8*rateMax;
//      numOfIters = max(numOfIters, 1);
//      numOfIters = 4;
//      double dtSrs = m_dt/numOfIters;
//      if (s_verbosity >= 3) {
//        pout() << "AMRLevelAdvectDiffuse::diffusiveAdvance: source " << m_level << endl;
//        pout() << "alpha = ";
//      }
//
//      double eps = 1e-10;
//      for (int isrs = 1; isrs <= numOfIters; isrs++) {
//        double alpha = (isrs-0.5)*dtSrs/m_dt;
//        if (s_verbosity >= 3)
//          pout() << alpha << " ";
//        for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//          FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//          FArrayBox rate(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//          FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//
//          ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
//          muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
//          getRate(rate, ETmp, muTemp, "ionization");
//
//          FArrayBox Uhs(m_UNew[dit()].box(), m_UNew.nComp()); // U at half step
//          Uhs.axby(m_UOld[dit()], m_UNew[dit()], 0.5, 0.5);
////          for (BoxIterator bit(m_grids.get(dit)); bit.ok(); ++bit) {
////            const IntVect& iv = bit();
////            if (abs(m_UOld[dit()](iv,0)) > eps && (abs(m_UNew[dit()](iv,0)-m_UOld[dit()](iv,0))/m_UOld[dit()](iv,0) > 0.5)) {
////              double gr = log(abs(m_UNew[dit()](iv,0)/m_UOld[dit()](iv,0)))/m_dt;
////              Uhs(iv, 0) = m_UOld[dit()](iv,0) * exp(gr * alpha*m_dt);
////            }
////            else
////              Uhs(iv, 0) = m_UOld[dit()](iv,0)*(1-alpha) + m_UNew[dit()](iv,0)*alpha;
////          }
//
////          Uhs.plus(dUTmp[dit()]);
//
//          m_srs[dit()].copy(rate);
//          m_srs[dit()].mult(Uhs);
//          m_srs[dit()].plus(m_phtzn.rate[dit()]);
//          Uhs.plus(m_srs[dit()], 0.5*m_dt);
//
//          m_srs[dit()].copy(rate);
//          m_srs[dit()].mult(Uhs);
//          m_srs[dit()].plus(m_phtzn.rate[dit()]);
//
//          dUTmp[dit()].plus(m_srs[dit()], m_dt * 1.0/numOfIters);
//        }
//      }
//
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        m_UNew[dit()].plus(dUTmp[dit()], m_grids[dit()], 0, 0);
//        m_ionNew[dit()].plus(dUTmp[dit()], m_grids[dit()], 0, 0);
//      }
//      if (s_verbosity >= 3)
//        pout() << endl;
//      break;
//    }
//
//    // Colella et al. 1999
//    case 5:
//    {
//      double alpha = 0.5;
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        m_srs[dit()].mult(alpha);
//        FArrayBox rate(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//        getRate(rate, m_field.m_Emag[dit()], m_mu[dit()], "ionization");
//        rate.mult(m_UNew[dit()]);
//        rate.plus(m_phtzn.rate[dit()]);
//        rate.mult(alpha);
//        m_srs[dit()].plus(rate);
//
//        // use this version so no update is made to ghost cells
//        m_UNew[dit()].plus(m_srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
//        m_ionNew[dit()].plus(m_srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
//      }
//
//      break;
//    }
      
    default:
      double alpha = 0.5;
      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateI(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateA(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        
        ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
        muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
        if (m_gas.m_uniformity) {
          getRate(rateI, ETmp, muTemp, "ionization");
          getRate(rateA, ETmp, muTemp, "attachment");
        } else {
          getRate(rateI, ETmp, m_neut[dit()], muTemp, "ionization");
          getRate(rateA, ETmp, m_neut[dit()], muTemp, "attachment");
        }
        
        srs[dit()].copy(rateI);
        srs[dit()].mult(m_UOld[dit()]);
        srs[dit()].plus(m_phtzn.rate[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        
        srs[dit()].copy(rateA);
        srs[dit()].mult(m_UOld[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], -m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 1);
      }
      
      break;
  }
  
  printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "AMRLevelAdvectDiffuse::diffusiveAdvanceEnd");
  printDiagnosticInfo (m_level, m_dx, m_grids, m_ionNew, "ion", "AMRLevelAdvectDiffuse::diffusiveAdvanceEnd");
  return newDt;
}

/*********/
Real
AMRLevelAdvectDiffuse::
strangAdvance(LevelData<FArrayBox>& a_diffusiveSrc)
{
  LevelFluxRegister* coarserFRPtr=NULL;
  LevelFluxRegister* finerFRPtr  =NULL;
  LevelData<FArrayBox>* coarserDataOldPtr = NULL;
  LevelData<FArrayBox>* coarserDataNewPtr = NULL;
  Real tCoarserOld, tCoarserNew;
  
  getCoarseDataPointers(&coarserDataOldPtr,
                        &coarserDataNewPtr,
                        &coarserFRPtr,
                        &finerFRPtr,
                        tCoarserOld, tCoarserNew);
  
  LevelData<FArrayBox> srs, srsTmp, srsRea, srsReaTmp; // divide sources into general source and reaction term
  srs.define(m_UNew);
  srsTmp.define(m_UNew);
  
  // call poissonSolve in case that the source charge at the current level is modified after last Poisson solve
  if (m_varyingField) {
    poissonSolve();
    fillMobility(false);
    fillAdvectionVelocity(false);
  }

  updateWithReactionContribution(m_UNew, m_field.m_Emag, m_mu, m_phtzn.rate, 0.5*m_dt);
  m_UNew.copyTo(m_UNew.interval(), m_UOld, m_UOld.interval());
  
  for (DataIterator dit=srs.dataIterator(); dit.ok(); ++dit) {
    srs[dit()].setVal(0.0);
    srsTmp[dit()].setVal(0.0);
  }

  for (DataIterator dit=srs.dataIterator(); dit.ok(); ++dit)
    srsTmp[dit()] += a_diffusiveSrc[dit()];
  srsTmp.exchange();
  
  //  patchgodunov L273
  //  localSource *= 0.5 * a_dt;
  // *finerFRPtr is set to zero at the beginning but *coarserFRPtr is incremented from prior value
  Real newDt = m_levelGodunov.step(m_UNew,
                                   m_flux,
                                   m_advVel,
                                   srsTmp,
                                   *coarserDataOldPtr,
                                   tCoarserOld,
                                   *coarserDataNewPtr,
                                   tCoarserNew,
                                   m_time,
                                   m_dt);
  
  if (m_hasFiner)
    finerFRPtr->setToZero();
  
  Interval interv(0, 0);
  if (m_hasDiffusion) {
    //compute Du = unew-uold and put uold back into unew
    //so we can advance using leveltga
    m_UNew.copyTo(interv, m_dU,   interv);
    for (DataIterator dit=m_dU.dataIterator(); dit.ok(); ++dit) {
      m_dU[dit()] -= m_UOld[dit()];
      m_dU[dit()] /= m_dt;
      //m_dU[dit()] += srs[dit()];
    }
    s_diffuseLevTGA->computeDiffusion(m_dUDiff, m_UOld, m_dU,
                                finerFRPtr, coarserFRPtr,
                                coarserDataOldPtr, coarserDataNewPtr,
                                m_time, tCoarserOld, tCoarserNew,
                                m_dt, m_level);
  }
  
  m_UOld.copyTo(interv, m_UNew, interv);
  
  if (m_doImplicitPoisson) {
    // m_flux is used in Poisson solve to set the variable coefficients
    if (m_varyingField)
      poissonSolveImplicit();
    for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
      Real eps = 1e-20;
      const Box& b = m_grids[dit()];
      FluxBox& curFlux = m_flux[dit()];
      FluxBox& advVel = m_advVel[dit()];
      advVel += eps;
      curFlux.divide(advVel, b, 0, 0);
    }
    fillAdvectionVelocity(true);
    for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
      const Box& b = m_grids[dit()];
      FluxBox& curFlux = m_flux[dit()];
      FluxBox& advVel = m_advVel[dit()];
      curFlux.mult(advVel, b, 0, 0);
    }
    // mobility is updated after the drift flux is computed so that the same drift flux is used in both Poisson's equation and the continuity equation
    fillMobility(true);
  }
    
  for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit)
    m_UNew[dit()].plus(m_dUDiff[dit()], m_dt);
  
  for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
    Interval UInterval(0, m_UNew.nComp()-1);
        
    Box curBox = m_grids.get(dit());
    FluxBox& curFlux = m_flux[dit()];
    m_dU[dit()].setVal(0.0);
    for (int idir = 0; idir < SpaceDim; idir++) {
      // Compute flux difference fHi-fLo
      FArrayBox diff(m_dU[dit()].box(), m_dU.nComp());
      diff.setVal(0.0);

      FORT_FLUXDIFFF(CHF_FRA(diff),
                     CHF_CONST_FRA(curFlux[idir]),
                     CHF_CONST_INT(idir),
                     CHF_BOX(curBox));

      // Add flux difference to dU
      m_dU[dit()] += diff;
    }
    m_dU[dit()] *= -m_dt/m_dx;
    m_UNew[dit()] += m_dU[dit()];
    
    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++) {
      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner)
        (*finerFRPtr).incrementCoarse(curFlux[idir],m_dt,dit(),
                                      UInterval,
                                      UInterval,idir);
      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser)
        (*coarserFRPtr).incrementFine(curFlux[idir],m_dt,dit(),
                                      UInterval,
                                      UInterval,idir);
    }
  }
  
  if (!m_doImplicitPoisson && m_varyingField)
    // Get the provisonal potential and mobility at the end of this step, which are used to fill ghost cells of the next finer level
    if (m_varyingField) {
      poissonSolve();
      fillMobility(true);
      fillAdvectionVelocity(true);
    }
  
  switch (m_sourceNumericalScheme) {
    
    // Dhali and Williams, 1987, but with the rate calculated using average field
    case 2: {
      double alpha = 0.5;
      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateI(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateA(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        
        ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
        muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
        if (m_gas.m_uniformity) {
          getRate(rateI, ETmp, muTemp, "ionization");
          getRate(rateA, ETmp, muTemp, "attachment");
        } else {
          getRate(rateI, ETmp, m_neut[dit()], muTemp, "ionization");
          getRate(rateA, ETmp, m_neut[dit()], muTemp, "attachment");
        }
        
        FArrayBox Uhs(m_UNew[dit()].box(), m_UNew.nComp()); // U at half step
        Uhs.axby(m_UOld[dit()], m_UNew[dit()], 1-alpha, alpha);
        srs[dit()].axby(rateI, rateA, 1.0, -1.0);
       
        srs[dit()].mult(m_UOld[dit()]);
        
        srs[dit()].plus(m_phtzn.rate[dit()]);
        Uhs.plus(srs[dit()], 0.5*m_dt);
        
        srs[dit()].copy(rateI);
        srs[dit()].mult(Uhs);
        srs[dit()].plus(m_phtzn.rate[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        
        srs[dit()].copy(rateA);
        srs[dit()].mult(Uhs);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], -m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 1);
      }
      
      break;
    }
     
    // differ with 2 only in Uhs, where source is calcuated using U mean
    case 3:
    {
      double alpha = 0.5;
      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateI(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateA(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        
        ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
        muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
        if (m_gas.m_uniformity) {
          getRate(rateI, ETmp, muTemp, "ionization");
          getRate(rateA, ETmp, muTemp, "attachment");
        } else {
          getRate(rateI, ETmp, m_neut[dit()], muTemp, "ionization");
          getRate(rateA, ETmp, m_neut[dit()], muTemp, "attachment");
        }
        
        FArrayBox Uhs(m_UNew[dit()].box(), m_UNew.nComp()); // U at half step
        Uhs.axby(m_UOld[dit()], m_UNew[dit()], 1-alpha, alpha);
        srs[dit()].axby(rateI, rateA, 1.0, -1.0);
        
        srs[dit()].mult(Uhs);
        
        srs[dit()].plus(m_phtzn.rate[dit()]);
        Uhs.plus(srs[dit()], 0.5*m_dt);
        
        srs[dit()].copy(rateI);
        srs[dit()].mult(Uhs);
        srs[dit()].plus(m_phtzn.rate[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        
        srs[dit()].copy(rateA);
        srs[dit()].mult(Uhs);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], -m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 1);
      }
      
      break;
    }
      
//    // LeVeque, p. 390, 2002
//    case 4: {
//      double rateMax = 0;
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        FArrayBox rateNew(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//        rateNew.setVal(0.0);
//        getRate(rateNew, m_field.m_Emag[dit()], m_mu[dit()], "ionization");
//        rateMax = max(rateMax, rateNew.max());
//      }
//
//      LevelData<FArrayBox> dUTmp;
//      dUTmp.define(m_UOld);
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        dUTmp[dit()].setVal(0.0);
//      }
//
//      int numOfIters = m_dt*8*rateMax;
//      numOfIters = max(numOfIters, 1);
//      numOfIters = 4;
//      double dtSrs = m_dt/numOfIters;
//      if (s_verbosity >= 3) {
//        pout() << "AMRLevelAdvectDiffuse::diffusiveAdvance: source " << m_level << endl;
//        pout() << "alpha = ";
//      }
//
//      double eps = 1e-10;
//      for (int isrs = 1; isrs <= numOfIters; isrs++) {
//        double alpha = (isrs-0.5)*dtSrs/m_dt;
//        if (s_verbosity >= 3)
//          pout() << alpha << " ";
//        for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//          FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//          FArrayBox rate(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//          FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//
//          ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
//          muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
//          getRate(rate, ETmp, muTemp, "ionization");
//
//          FArrayBox Uhs(m_UNew[dit()].box(), m_UNew.nComp()); // U at half step
//          Uhs.axby(m_UOld[dit()], m_UNew[dit()], 0.5, 0.5);
////          for (BoxIterator bit(m_grids.get(dit)); bit.ok(); ++bit) {
////            const IntVect& iv = bit();
////            if (abs(m_UOld[dit()](iv,0)) > eps && (abs(m_UNew[dit()](iv,0)-m_UOld[dit()](iv,0))/m_UOld[dit()](iv,0) > 0.5)) {
////              double gr = log(abs(m_UNew[dit()](iv,0)/m_UOld[dit()](iv,0)))/m_dt;
////              Uhs(iv, 0) = m_UOld[dit()](iv,0) * exp(gr * alpha*m_dt);
////            }
////            else
////              Uhs(iv, 0) = m_UOld[dit()](iv,0)*(1-alpha) + m_UNew[dit()](iv,0)*alpha;
////          }
//
////          Uhs.plus(dUTmp[dit()]);
//
//          m_srs[dit()].copy(rate);
//          m_srs[dit()].mult(Uhs);
//          m_srs[dit()].plus(m_phtzn.rate[dit()]);
//          Uhs.plus(m_srs[dit()], 0.5*m_dt);
//
//          m_srs[dit()].copy(rate);
//          m_srs[dit()].mult(Uhs);
//          m_srs[dit()].plus(m_phtzn.rate[dit()]);
//
//          dUTmp[dit()].plus(m_srs[dit()], m_dt * 1.0/numOfIters);
//        }
//      }
//
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        m_UNew[dit()].plus(dUTmp[dit()], m_grids[dit()], 0, 0);
//        m_ionNew[dit()].plus(dUTmp[dit()], m_grids[dit()], 0, 0);
//      }
//      if (s_verbosity >= 3)
//        pout() << endl;
//      break;
//    }
//
//    // Colella et al. 1999
//    case 5:
//    {
//      double alpha = 0.5;
//      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
//        m_srs[dit()].mult(alpha);
//        FArrayBox rate(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
//        getRate(rate, m_field.m_Emag[dit()], m_mu[dit()], "ionization");
//        rate.mult(m_UNew[dit()]);
//        rate.plus(m_phtzn.rate[dit()]);
//        rate.mult(alpha);
//        m_srs[dit()].plus(rate);
//
//        // use this version so no update is made to ghost cells
//        m_UNew[dit()].plus(m_srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
//        m_ionNew[dit()].plus(m_srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
//      }
//
//      break;
//    }
      
    default:
      double alpha = 0.5;
      for (DataIterator dit=m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox ETmp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateI(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox rateA(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        FArrayBox muTemp(m_field.m_Emag[dit()].box(), m_field.m_Emag.nComp());
        
        ETmp.axby(m_fieldOld.m_Emag[dit()], m_field.m_Emag[dit()], 1-alpha, alpha);
        muTemp.axby(m_muOld[dit()], m_mu[dit()], 1-alpha, alpha);
        if (m_gas.m_uniformity) {
          getRate(rateI, ETmp, muTemp, "ionization");
          getRate(rateA, ETmp, muTemp, "attachment");
        } else {
          getRate(rateI, ETmp, m_neut[dit()], muTemp, "ionization");
          getRate(rateA, ETmp, m_neut[dit()], muTemp, "attachment");
        }
        
        srs[dit()].copy(rateI);
        srs[dit()].mult(m_UOld[dit()]);
        srs[dit()].plus(m_phtzn.rate[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 0);
        
        srs[dit()].copy(rateA);
        srs[dit()].mult(m_UOld[dit()]);
        // use this version so no update is made to ghost cells
        m_UNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], -m_dt, 0, 0);
        m_ionNew[dit()].plus(srs[dit()], m_grids[dit()], m_grids[dit()], m_dt, 0, 1);
      }
      
      break;
  }
  
  printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "AMRLevelAdvectDiffuse::diffusiveAdvanceEnd");
  printDiagnosticInfo (m_level, m_dx, m_grids, m_ionNew, "ion", "AMRLevelAdvectDiffuse::diffusiveAdvanceEnd");
  return newDt;
}

/*********/
void
AMRLevelAdvectDiffuse::
getRate(LevelData<FArrayBox>& a_rate, const LevelData<FArrayBox>& a_Emag, const LevelData<FArrayBox>& a_mu, string processName) {
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    getRate(a_rate[dit()], a_Emag[dit()], a_mu[dit()], processName);
  }
  a_rate.exchange();
}

/*********/
void
AMRLevelAdvectDiffuse::
getRate(FArrayBox& a_rate, const FArrayBox& a_Emag, const FArrayBox& a_mu, string processName) {
  a_rate.setVal(0);
  a_rate.copy(a_mu);
  // store ve
  a_rate.mult(a_Emag);
  for (BoxIterator bit(a_rate.box()); bit.ok(); ++bit) {
    const IntVect& iv = bit();
    Real E, a;
    Real n = m_gas.m_N;
    E = a_Emag(iv, 0);
    a = m_gas.EDpdentProcs[processName].value(E/n)*n;
    a_rate(iv, 0) *= a;
  }
  
}

/*********/
void
AMRLevelAdvectDiffuse::
getRate(LevelData<FArrayBox>& a_rate, const LevelData<FArrayBox>& a_Emag, const LevelData<FArrayBox>& a_n, const LevelData<FArrayBox>& a_mu, string processName) {
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    getRate(a_rate[dit()], a_Emag[dit()], a_n[dit()], a_mu[dit()], processName);
  }
  a_rate.exchange();
}

/*********/
void
AMRLevelAdvectDiffuse::
getRate(FArrayBox& a_rate, const FArrayBox& a_Emag, const FArrayBox& a_n, const FArrayBox& a_mu, string processName) {
  a_rate.setVal(0);
  a_rate.copy(a_mu);
  // store ve
  a_rate.mult(a_Emag);
  for (BoxIterator bit(a_rate.box()); bit.ok(); ++bit) {
    const IntVect& iv = bit();
    Real E, a;
    Real n = a_n(iv, 0);
    E = a_Emag(iv, 0);
    a = m_gas.EDpdentProcs[processName].value(E/n)*n;
    a_rate(iv, 0) *= a;
  }
  
}

/*******/
void
AMRLevelAdvectDiffuse::
setSolverCoef(Real a_alpha, Real a_beta)
{
  // now set alpha and beta of all the operators in the solver.
  // this includes resetting the relaxation coefficient
  Vector<MGLevelOp<LevelData<FArrayBox> >* > ops = s_diffuseAMRMG->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
  {
    LevelTGAHelmOp<FArrayBox,FluxBox>* helmop = (LevelTGAHelmOp<FArrayBox,FluxBox>*) ops[iop];
    helmop->setAlphaAndBeta(a_alpha, a_beta);
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
postTimeStep()
{
  if (s_verbosity >= 2) {
    pout() << "AMRLevelAdvectDiffuse::postTimeStep " << m_level << endl;
  }
  if (m_hasFiner) {
    if (m_doImplicitReflux) {
      doImplicitReflux();
      if (m_varyingField)
        poissonSolveComposite();
    }
    // explicit Reflux
    else {
      if (m_doImplicitPoisson) {
        if (m_varyingField)
          poissonSolveImplicitComposite();
        
      } else {
        Real scale = -1.0/m_dx;
        m_fluxRegister.reflux(m_UNew,scale);
        // Average from finer level data
        AMRLevelAdvectDiffuse* amrGodFinerPtr = getFinerLevel();
        if (m_varyingField)
          poissonSolveComposite();
      }
    } // end if we're doing explicit refluxing
  } // end if there is a finer level

  Real time_eps = 1.0e-10;
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps)) {
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    int finest_level = hierarchy.size()-1;
    
    for (int lev = m_level; lev <= finest_level; lev++) {
      if (lev == m_level) {
        hierarchy[lev]->fillMobility(true);
        hierarchy[lev]->fillAdvectionVelocity(true);
      }
      else {
        hierarchy[lev]->fillMobility(false);
        hierarchy[lev]->fillAdvectionVelocity(false);
      }
      
      
//      LevelData<FArrayBox> srsOld;
//      srsOld.define(hierarchy[lev]->m_srs);
//      for (DataIterator dit=grids[lev].dataIterator(); dit.ok(); ++dit)
//        srsOld[dit()].setVal(0.0);
//      hierarchy[lev]->m_srs.copyTo(srsOld);
//      hierarchy[lev]->makeSource(hierarchy[lev]->m_UNew, 0.0);
//      for (DataIterator dit=grids[lev].dataIterator(); dit.ok(); ++dit) {
//        hierarchy[lev]->m_UNew[dit()].plus(srsOld[dit()], -0.5*hierarchy[lev]->m_dt);
//        hierarchy[lev]->m_ionNew[dit()].plus(srsOld[dit()], -0.5*hierarchy[lev]->m_dt);
//        hierarchy[lev]->m_UNew[dit()].plus(hierarchy[lev]->m_srs[dit()], 0.5*hierarchy[lev]->m_dt);
//        hierarchy[lev]->m_ionNew[dit()].plus(hierarchy[lev]->m_srs[dit()], 0.5*hierarchy[lev]->m_dt);
//      }
    }
    
    // Average from finer level data
    for (int lev = finest_level; lev >= m_level; lev--)
      if (hierarchy[lev]->m_hasCoarser) {
        hierarchy[lev]->m_coarseAverage.averageToCoarse(hierarchy[lev-1]->m_UNew, hierarchy[lev]->m_UNew);
        hierarchy[lev]->m_coarseAverageM.averageToCoarse(hierarchy[lev-1]->m_ionNew, hierarchy[lev]->m_ionNew);
      }
  }
  
  if (m_phtzn.runSolve)
    photoionizationSolve();
  
  outputStepStats(AMPLIFIOut);
  printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "AMRLevelAdvectDiffuse::postTimeStep");
  printDiagnosticInfo (m_level, m_dx, m_grids, m_ionNew, "ion", "AMRLevelAdvectDiffuse::postTimeStep");
  //  outputDataForCheck (m_level, m_grids, m_field.m_E);  
}

/*******/
void
AMRLevelAdvectDiffuse::
doImplicitReflux()
{
  // first find out what time coarser level is (if it exists)
  Real time_eps = 1.0e-10;
  // do multilevel operations if this is the coarsest level or if
  // coarser level is not at same time as this level. otherwise,
  // defer this until we get down to the coarsest level at this time.
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps))
  {
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    
    int finest_level = hierarchy.size()-1;
    
    // now do implicit refluxing
    // Vector of pointers to LevelData of FABS
    Vector<LevelData<FArrayBox>* > refluxCor(finest_level+1, NULL);
    Vector<LevelData<FArrayBox>* > refluxRHS(finest_level+1, NULL);
    // collectUN: AMR vector containing soln at new time
    Vector<LevelData<FArrayBox>* > collectUN(finest_level+1, NULL);
    
    // loop over levels, allocate storage, set up for AMRSolve
    // if coarser level exists, define it as well for BCs.
    int startLev = Max(m_level-1, 0);
    
    for (int lev = startLev; lev<= finest_level; lev++)
    {
      // rhs has no ghost cells, correction does
      refluxRHS[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      refluxCor[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Unit);
      collectUN[lev]  = &(hierarchy[lev]->m_UNew);
      for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
      {
        (*refluxRHS[lev])[dit()].setVal(0.0);
        (*refluxCor[lev])[dit()].setVal(0.0);
      }
    } // end loop over levels for setup.
    
    // now loop over levels and set up RHS
    // note that we don't need to look at finest level here,
    // since it will have no finer level to get a reflux correction
    // from.   Also this starts at m_level instead of startLev since,
    // in the case m_level > 0, m_level-1 is only used for boundary conditions
    for (int lev= m_level; lev < finest_level; lev++)
    {
      // see line 692 of BaseLevelTGA.H
      // Now increment the flux registers -- note that because of the way
      // we defined the fluxes, the dt multiplier is already in the
      // flux.
      Real dxLev = hierarchy[lev]->m_dx;
      Real refluxScale = 1.0/dxLev;
      
      hierarchy[lev]->m_fluxRegister.reflux(*refluxRHS[lev], refluxScale);
    }
    
    
    int lbase = m_level;
    int lmax  = finest_level;
    // this resets the coeffients including eta, alpha, beta
    Real alpha = 1.0;
    Real beta  = -m_dt;
    setSolverCoef(alpha, beta);
    s_diffuseAMRMG->solve(refluxCor, refluxRHS, lmax, lbase);
    
    for (int lev= m_level; lev <= finest_level; lev++)
    {
      for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
      {
        hierarchy[lev]->m_UNew[dit()] += (*refluxCor[lev])[dit()];
      }
    }
    
    //remember that startLev can be different from m_level
    for (int lev = startLev; lev<= finest_level; lev++)
    {
      delete refluxRHS[lev];
      delete refluxCor[lev];
      refluxRHS[lev] = NULL;
      refluxCor[lev] = NULL;
    }
  } //end if times are in alignment
}

/*******/
void
AMRLevelAdvectDiffuse::
poissonSolve() {

  pout() << "AMRLevelAdvectDiffuse::poissonSolve " << m_level << endl;
    
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    
    int finest_level = hierarchy.size()-1;
    
    // solve Poisson's equation
    // Vector of pointers to LevelData of FABS
    Vector<LevelData<FArrayBox>* > phi(finest_level+1, NULL);
    Vector<LevelData<FArrayBox>* > rhs(finest_level+1, NULL);
    
    // loop over levels, allocate storage, set up for AMRSolve
    // if coarser level exists, define it as well for BCs.
    int startLev = Max(m_level-1, 0);
    
    RefCountedPtr<AMRPoissonOp> amrpop = RefCountedPtr<AMRPoissonOp>((AMRPoissonOp*)s_EPotOpFact->AMRnewOp(m_problem_domain));
    
    for (int lev = startLev; lev<= m_level; lev++) {
      // rhs has no ghost cells, phi does
      rhs[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      phi[lev]  = new LevelData<FArrayBox>(grids[lev], 1, m_numGhost*IntVect::Unit);
      
      hierarchy[lev]->m_ionNew.copyTo(Interval(0, 0), *rhs[lev], Interval(0, 0));
      for (DataIterator dit=rhs[lev]->dataIterator(); dit.ok(); ++dit) {
        (*rhs[lev])[dit()].minus(hierarchy[lev]->m_ionNew[dit()], 1, 0);
        (*rhs[lev])[dit()] -= hierarchy[lev]->m_UNew[dit()];
      }
      
      amrpop->scale(*rhs[lev], -1);
            
      hierarchy[lev]->m_phi.copyTo(*phi[lev]);
    } // end loop over levels for setup.
    
    if (m_hasCoarser) {
      AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
//      interpolateInTime(*phi[startLev], amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_phi, m_time+0.5*m_dt, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_time);
      interpolateInTime(*phi[startLev], amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_phi, m_time+m_dt, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_time);
    }
    
    int lbase = m_level;
    int lmax  = m_level;
    s_EPotAMRMG->solve(phi, rhs, lmax, lbase, false);
        
    (*phi[m_level]).copyTo(m_phi);
    
    computeEField(true);
  
    //remember that startLev can be different from m_level
    for (int lev = startLev; lev<= finest_level; lev++) {
      delete rhs[lev];
      delete phi[lev];
      rhs[lev] = NULL;
      phi[lev] = NULL;
    }
}


/*******/
void
AMRLevelAdvectDiffuse::
poissonSolveComposite() {
  // first find out what time coarser level is (if it exists)
  Real time_eps = 1.0e-10;
  // do multilevel operations if this is the coarsest level or if
  // coarser level is not at same time as this level. otherwise,
  // defer this until we get down to the coarsest level at this time.
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps)) {
    if (s_verbosity >= 3)
      pout() << "AMRLevelAdvectDiffuse::poissonSolveComposite " << m_level << endl;
    
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    
    int finest_level = hierarchy.size()-1;
    
    // solve Poisson's equation
    // Vector of pointers to LevelData of FABS
    Vector<LevelData<FArrayBox>* > phi(finest_level+1, NULL);
    Vector<LevelData<FArrayBox>* > rhs(finest_level+1, NULL);
    
    // loop over levels, allocate storage, set up for AMRSolve
    // if coarser level exists, define it as well for BCs.
    int startLev = Max(m_level-1, 0);
    
    RefCountedPtr<AMRPoissonOp> amrpop = RefCountedPtr<AMRPoissonOp>((AMRPoissonOp*)s_EPotOpFact->AMRnewOp(m_problem_domain));
    
    for (int lev = startLev; lev<= finest_level; lev++) {
      // rhs has no ghost cells, phi does
      rhs[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      phi[lev]  = new LevelData<FArrayBox>(grids[lev], 1, m_numGhost*IntVect::Unit);
    
      hierarchy[lev]->m_ionNew.copyTo(Interval(0, 0), *rhs[lev], Interval(0, 0));
      for (DataIterator dit=rhs[lev]->dataIterator(); dit.ok(); ++dit) {
        (*rhs[lev])[dit()].minus(hierarchy[lev]->m_ionNew[dit()], 1, 0);
        (*rhs[lev])[dit()] -= hierarchy[lev]->m_UNew[dit()];
      }
      
      amrpop->scale(*rhs[lev], -1);
      
      hierarchy[lev]->m_phi.copyTo(*phi[lev]);
      
    } // end loop over levels for setup.
    
    if (m_hasCoarser) {
      AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
//      interpolateInTime(*phi[startLev], amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_phi, m_time+0.5*m_dt, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_time);
      interpolateInTime(*phi[startLev], amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_phi, m_time+m_dt, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_time);
    }
    
    int lbase = m_level;
    int lmax  = finest_level;
    s_EPotAMRMG->solve(phi, rhs, lmax, lbase, false);
    
    for (int lev = m_level; lev <= finest_level; lev++)
      (*phi[lev]).copyTo(hierarchy[lev]->m_phi);
    
    // need to find fields for m_level to finest_level, because they are required to calculate the photoionziation source
    
    // multigrid smoothes down first and then up; fine level solution is more accurate! This is really necessary; otherwise, there would be discontinuity at the coarse-fine interface!
    for (int lev = finest_level; lev >= m_level; lev--)
      if (hierarchy[lev]->m_hasCoarser)
        hierarchy[lev]->m_coarseAverage.averageToCoarse(hierarchy[lev-1]->m_phi, hierarchy[lev]->m_phi);
    
    computeEField(true);
    for (int lev = m_level+1; lev <= finest_level; lev++)
      hierarchy[lev]->computeEField(false);

    //remember that startLev can be different from m_level
    for (int lev = startLev; lev<= finest_level; lev++) {
      delete rhs[lev];
      delete phi[lev];
      rhs[lev] = NULL;
      phi[lev] = NULL;
    }
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
poissonSolveImplicit() {
  
  pout() << "AMRLevelAdvectDiffuse::poissonSolveImplicit " << m_level << endl;
  
  Vector<AMRLevelAdvectDiffuse*>         hierarchy;
  Vector<int>                            refRat;
  Vector<DisjointBoxLayout>              grids;
  Real                                   lev0Dx;
  ProblemDomain                          lev0Domain;
  getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
  
  int finest_level = hierarchy.size()-1;
  
  // solve Poisson's equation
  // Vector of pointers to LevelData of FABS
  Vector<LevelData<FArrayBox>* > phi(finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > rhs(finest_level+1, NULL);

  // loop over levels, allocate storage, set up for AMRSolve
  // if coarser level exists, define it as well for BCs.
  int startLev = Max(m_level-1, 0);
  
  for (int lev = startLev; lev<= m_level; lev++) {
    // rhs has no ghost cells, phi does
    rhs[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
    phi[lev]  = new LevelData<FArrayBox>(grids[lev], 1, m_numGhost*IntVect::Unit);
    hierarchy[lev]->m_phi.copyTo(*phi[lev]);
  }
  // end loop over levels for setup.
  
  setPoissonCoeffAB();
  
  int lev = m_level;
  AMRLevelAdvectDiffuse* cl = hierarchy[lev];
  for (DataIterator dit=rhs[lev]->dataIterator(); dit.ok(); ++dit) {
    (*rhs[lev])[dit()].setVal(0.0);
  }
  
  if(m_hasFiner) {
    AMRLevelAdvectDiffuse* fl = hierarchy[lev+1];
    // the first 2 terms in the RHS of eq. (4.49) in Colella et al. (1999b) give bCoeff, which contains dt_c.
    cl->getDivDeltaFlux((*rhs[lev]), &(cl->m_EEdgeOld), &(fl->m_EEdgeCoarseOld), s_bCoef[lev], 1.0/cl->m_dx);
  }
  for (DataIterator dit=grids[lev].dataIterator(); dit.ok(); ++dit) {
    (*rhs[lev])[dit()] += cl->m_UNew[dit()];
    (*rhs[lev])[dit()].minus(hierarchy[lev]->m_ionNew[dit()], 0, 0);
    (*rhs[lev])[dit()].plus(hierarchy[lev]->m_ionNew[dit()], 1, 0);
    
    //see eq. (3.63) in Colella et al. (1999a) or (3.23) in Colella et al. (1999b). Gamma_diff contains minus sign.
    (*rhs[lev])[dit()].plus(cl->m_dUDiff[dit()], m_dt);
  }
  
  if (m_hasCoarser) {
    AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
//      interpolateInTime(*phi[startLev], amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_phi, m_time+0.5*m_dt, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_time);
    interpolateInTime(*phi[startLev], amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_phi, m_time+m_dt, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_time);
  }
  
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
  amrpopEPotVC->define(lev0Domain, grids, refRat, lev0Dx, m_EPotbcFunc, 0.0, s_aCoef, -1.0, s_bCoef);
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > EPotImpOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(amrpopEPotVC);
  RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > EPotImpAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >
  (new AMRMultiGrid<LevelData<FArrayBox> >());
  EPotImpAMRMG->define(lev0Domain, *EPotImpOpFact, &s_EPotImpBotSolver, hierarchy.size());
  EPotImpAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
                                      maxIter, tolerance, hang, normThresh);
  EPotImpAMRMG->m_verbosity = mgverb;
  s_EPotImpBotSolver.m_verbosity = mgverb;
    
  int lbase = m_level;
  int lmax  = m_level;
  
  EPotImpAMRMG->solve(phi, rhs, lmax, lbase, false);
  
  (*phi[m_level]).copyTo(m_phi);
  
  computeEField(true);
  
  //remember that startLev can be different from m_level
  for (int lev = startLev; lev<= finest_level; lev++) {
    delete rhs[lev];
    delete phi[lev];
    rhs[lev] = NULL;
    phi[lev] = NULL;
  }
}

void
AMRLevelAdvectDiffuse::
setPoissonCoeffAB() {
  LevelData<FluxBox>&  muEdge = m_muEdge;
  
  Real eps = 1e-20;
  m_flux.copyTo(*s_bCoef[m_level]);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    (*s_aCoef[m_level])[dit()].setVal(0.0);
    
    FluxBox advVel;
    advVel.define(m_advVel[dit()]);
    advVel.copy(m_advVel[dit()]);
    advVel += eps;
    const Box& b = m_grids[dit()];
    (*s_bCoef[m_level])[dit()].divide(advVel, b, 0, 0);
    (*s_bCoef[m_level])[dit()].mult(muEdge[dit()], b, 0, 0);
    (*s_bCoef[m_level])[dit()] *= m_dt;
    (*s_bCoef[m_level])[dit()] += 1.0;
  }
  (*s_aCoef[m_level]).exchange();
  (*s_bCoef[m_level]).exchange();
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::setPoissonCoeffAB " << m_level << endl;
  if (s_verbosity >= 3) {
     printDiagnosticInfo (m_level, m_dx, m_grids, *s_bCoef[m_level], "bCoeff", "setPoissonCoeffAB");
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
poissonSolveImplicitComposite() {
  // first find out what time coarser level is (if it exists)
  Real time_eps = 1.0e-10;
  // do multilevel operations if this is the coarsest level or if
  // coarser level is not at same time as this level. otherwise,
  // defer this until we get down to the coarsest level at this time.
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps)) {
    if (s_verbosity >= 3)
      pout() << "AMRLevelAdvectDiffuse::poissonSolveImplicitComposite " << m_level << endl;
      
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
      hierarchy[lev]->setPoissonCoeffABComposite(m_dt);
    
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
    
    for (int lev = finest_level; lev>= m_level; lev--) {
      if (hierarchy[lev]->m_hasFiner && lev >= m_level) {
        AMRLevelAdvectDiffuse* cl = hierarchy[lev];
        AMRLevelAdvectDiffuse* fl = hierarchy[lev+1];
        
        // revisit: For cases with >2 levels, is dt here the one for lev or lbase?
        // Define coe = dt*neMidStep*mu
        LevelData<FluxBox> coe;
//        Real scale = m_dt/cl->m_dt;
        Real scale = 1.0;
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
    if (s_verbosity >= 3)
      pout() << "iter = " << iter << " " << "totVolCharge = " << totVolCharge << "  " << "totSurfCharge = " << totSurfCharge << endl;
    
    // save the predictor solution
    for (int lev = m_level; lev <= finest_level; lev++)
      hierarchy[lev]->m_phi.copyTo(hierarchy[lev]->m_phiOld);
    
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
        
        // need to find fields for m_level to finest_level, because they are required to calculate the photoionziation source

        // multigrid smoothes down first and then up; fine level solution is more accurate! This is really necessary; otherwise, there would be discontinuity at the coarse-fine interface!
        for (int lev = finest_level; lev >= m_level; lev--)
          if (hierarchy[lev]->m_hasCoarser)
            hierarchy[lev]->m_coarseAverage.averageToCoarse(hierarchy[lev-1]->m_phi, hierarchy[lev]->m_phi);
        
        computeEField(true);
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
setPoissonCoeffABComposite(Real commonDt) {
  
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
    
  if (m_hasFiner) {
    AMRLevelAdvectDiffuse* amrGodFinerPtr = getFinerLevel();
//    (*s_bCoefComp[m_level+1]).exchange();
    amrGodFinerPtr->m_coarseAverageFace.averageToCoarse(*s_bCoefComp[m_level], *s_bCoefComp[m_level+1]);

  }
  
  if (s_verbosity >= 3) {
     printDiagnosticInfo (m_level, m_dx, m_grids, *s_bCoefComp[m_level], "bCoeffComp", "setPoissonCoeffABComp");
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
getDivDeltaFlux(LevelData<FArrayBox>& rhs, LevelData<FluxBox>* edgeField, LevelData<FluxBox>* edgeFieldFine, LevelData<FluxBox>* coeff, Real scale) {
  AMRLevelAdvectDiffuse* fl = getFinerLevel();
  LevelFluxRegister frDEFlux;
  frDEFlux.define(fl->m_grids, m_grids, fl->m_problem_domain, m_ref_ratio, 1);
  frDEFlux.setToZero();
  
  AMRLevelAdvectDiffuse* l = this;
  for (DataIterator dit=(*l).m_grids.dataIterator(); dit.ok(); ++dit) {
    Interval UInterval(0, 0);
    FluxBox& curEEdge = (*edgeField)[dit()];
    FluxBox dEFluxCoarse;
    dEFluxCoarse.define(curEEdge);
    dEFluxCoarse.copy(curEEdge);
    if (coeff != NULL) {
      FluxBox& curCoef = (*coeff)[dit()];
      // curCoef has no ghost cells
      dEFluxCoarse.mult(curCoef, curCoef.box(), 0, 0);
    }
    for (int idir = 0; idir < SpaceDim; idir++) {
      frDEFlux.incrementCoarse(dEFluxCoarse[idir], 1.0, dit(), UInterval, UInterval, idir);
    }
  }
  
  LevelData<FluxBox> coefInterpToFine;
  if (coeff != NULL) {
    coefInterpToFine.define(*edgeFieldFine);
    fl->m_fineInterpFace.interpToFine(coefInterpToFine, *coeff);
  }
  l = fl;
  for (DataIterator dit=(*l).m_grids.dataIterator(); dit.ok(); ++dit) {
    Interval UInterval(0, 0);
    FluxBox& curEEdge = (*edgeFieldFine)[dit()];
    FluxBox dEFluxFine;
    dEFluxFine.define(curEEdge);
    dEFluxFine.copy(curEEdge);
    if (coeff != NULL) {
      FluxBox& curCoef = coefInterpToFine[dit()];
      dEFluxFine.mult(curCoef, curCoef.box(), 0, 0);
    }
    for (int idir = 0; idir < SpaceDim; idir++) {
      frDEFlux.incrementFine(dEFluxFine[idir], 1.0, dit(), UInterval, UInterval, idir);
    }
  }
  
  frDEFlux.reflux(rhs, scale);
}

void AMRLevelAdvectDiffuse::
computeEField(bool timeInterpForGhost)
{
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::computeEField " << m_level << endl;
  
 
  if (m_hasCoarser) {
    AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
//    Real eps = 1.0e-10;
//    if (timeInterpForGhost)
//      fillGhost(m_pwl, m_phi, amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_phi, amrGodCoarserPtr->m_time, m_dt, m_time+m_dt);
//    else
//      fillGhost(m_pwl, m_phi, amrGodCoarserPtr->m_phi, 0, amrGodCoarserPtr->m_phi, 0, eps, 0);
    // need to use quadratic interpolation; otherwise, a fluctuation appears at the coarse-fine interface.
    if (timeInterpForGhost)
      fillGhostQuad(m_qcfi, m_phi, amrGodCoarserPtr->m_phiOld, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_phi, amrGodCoarserPtr->m_time, m_time+m_dt);
    else
      m_qcfi.coarseFineInterp(m_phi, amrGodCoarserPtr->m_phi);
  }
  

  // Set the ghost cells outside the physical boundary
  // Set by assuming the same potential difference as the neighbors
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    const Box& b = m_grids[dit()];
    FArrayBox& phiFab = m_phi[dit()];
    
    if (!m_problem_domain.domainBox().contains(phiFab.box())) {
      for (int dir = 0; dir < SpaceDim; ++dir) {
        const Box bLoGhost = adjCellBox(b, dir, Side::Lo, 1);
        const Box bHiGhost = adjCellBox(b, dir, Side::Hi, 1);
        Box bCenter1 = b & grow(m_problem_domain,-BASISV(dir));
        Box bCenter2 = b & grow(m_problem_domain,-BASISV(dir)*2);
        if (!m_problem_domain.domainBox().contains(bLoGhost)) {
          const Box bLo1     = b & adjCellLo(bCenter1,dir);
          const Box bLo2     = b & adjCellLo(bCenter2,dir);
          phiFab.copy(phiFab, bLo1, 0, bLoGhost, 0, m_phi.nComp());
          phiFab.mult(2.0, bLoGhost);
          phiFab.minus(phiFab, bLo2, bLoGhost, 0, 0);
        }
        if (!m_problem_domain.domainBox().contains(bHiGhost)) {
          const Box bHi1     = b & adjCellHi(bCenter1,dir);
          const Box bHi2     = b & adjCellHi(bCenter2,dir);
          phiFab.copy(phiFab, bHi1, 0, bHiGhost, 0, m_phi.nComp());
          phiFab.mult(2.0, bHiGhost);
          phiFab.minus(phiFab, bHi2, bHiGhost, 0, 0);
        }
      }
    }
  }
  m_phi.exchange();

    
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    m_field.m_EEdge[dit()].setVal(0.0);
  }
  
  Gradient::levelGradientMAC(m_field.m_EEdge, m_phi, m_dx);
  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    m_field.m_EEdge[dit()].negate();
  }
  
  EdgeToCell(m_field.m_EEdge, m_field.m_E);
  
  // find magnitude
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    FArrayBox& gradFab = m_field.m_E[dit()];
    FArrayBox& gradMagFab = m_field.m_Emag[dit()];
    const Box& b = gradMagFab.box();
    FORT_MAGNITUDEF(CHF_FRA1(gradMagFab,0),
                    CHF_CONST_FRA(gradFab),
                    CHF_BOX(b));
  }

  if (s_verbosity >= 3) {
    printDiagnosticInfo (m_level, m_dx, m_grids, m_field.m_E, "E", "computeEField");
    printDiagnosticInfo (m_level, m_dx, m_grids, m_field.m_Emag, "Emag", "computeEField");
    printDiagnosticInfo (m_level, m_dx, m_grids, m_field.m_EEdge, "EEdge", "computeEField");
  }

}

  
/*******/
void
AMRLevelAdvectDiffuse::
photoionizationSolve() {
  // first find out what time coarser level is (if it exists)
  Real time_eps = 1.0e-10;
  // do multilevel operations if this is the coarsest level or if
  // coarser level is not at same time as this level. otherwise,
  // defer this until we get down to the coarsest level at this time.
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps)) {
    if (s_verbosity >= 3)
      pout() << "AMRLevelAdvectDiffuse::photoionizationSolve " << m_level << endl;
    
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    
    int finest_level = hierarchy.size()-1;
    
    // solve Poisson's equation
    // Vector of pointers to LevelData of FABS
    Vector<LevelData<FArrayBox>* > phi(finest_level+1, NULL);
    Vector<LevelData<FArrayBox>* > rhs(finest_level+1, NULL);
    
    // loop over levels, allocate storage, set up for AMRSolve
    // if coarser level exists, define it as well for BCs.
    int startLev = Max(m_level-1, 0);
    
    Vector<MGLevelOp<LevelData<FArrayBox> >* > ops = s_PIAMRMG->getAllOperators();
    
    for (int lev = startLev; lev<= finest_level; lev++) {
      // rhs has no ghost cells, while phi does.
      rhs[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      phi[lev]  = new LevelData<FArrayBox>(grids[lev], 1, m_numGhost*IntVect::Unit);
      
      hierarchy[lev]->m_mu.copyTo(*(rhs[lev]));
      for (DataIterator dit =  hierarchy[lev]->m_grids.dataIterator(); dit.ok(); ++dit) {
        FArrayBox& rhsFab = (*rhs[lev])[dit()];
        FArrayBox& EmagFab = hierarchy[lev]->m_field.m_Emag[dit()];
        rhsFab.mult(EmagFab);
        rhsFab.mult(hierarchy[lev]->m_UNew[dit()]);
        rhsFab.mult(quenchingFactor*PIXi);
        for (BoxIterator bit(hierarchy[lev]->m_grids.get(dit)); bit.ok(); ++bit) {
          const IntVect& iv = bit();
          Real E, ai;
          Real n;
          E = EmagFab(iv, 0);
          if (m_gas.m_uniformity)
            n = m_gas.m_N;
          else
            n = (hierarchy[lev]->m_neut[dit()])(iv,0);
          
          ai = m_gas.EDpdentProcs["ionization"].value(E/n)*n;
          rhsFab(iv, 0) *= ai;
        }
      }
      
    }
    
    ParmParse ppSolver("PISolver");
    bool PITesting;
    ppSolver.query("testing", PITesting);
    if (PITesting)
      PISetTestRHS(rhs, refRat, lev0Domain, lev0Dx, startLev, finest_level, PITestProbType::gaussian);
    
//    PISetTestRHS(rhs, refRat, lev0Domain, lev0Dx, startLev, finest_level, PITestProbType::gaussian);
    
//    for (int lev = startLev; lev<= finest_level; lev++)
//      if (s_verbosity >= 3)
//        printDiagnosticInfo (lev, hierarchy[lev]->m_dx, hierarchy[lev]->m_grids, *rhs[lev], "rhs", "AMRLevelAdvectDiffuse::photoionizationSolve");
    
    for (int i = 0; i < m_phtzn.ncomps; i++) {
      for (int iop = 0; iop < ops.size(); iop++) {
        AMRPoissonOp* amrpop = (AMRPoissonOp*) ops[iop];
        amrpop->setAlphaAndBeta(-3*m_phtzn.lambda[i]*m_phtzn.lambda[i], 1.0);
      }
      
      for (int lev = startLev; lev<= finest_level; lev++) {
        // need to do time interpolation for m_level-1??? no, Poisson solve at a finer level is done earlier than the next coarsers level. There is no point to do time interpolation.
        // no need for rhs because of no ghost cells
        hierarchy[lev]->m_phtzn.Psi.copyTo(Interval(i,i), *phi[lev], Interval(0,0));
        
      } // end loop over levels for setup.
      
      if (i == 0) {
        for (int lev = startLev; lev<= finest_level; lev++)
          (s_PIAMRMG->levelOp(lev)).scale(*rhs[lev], -3*m_phtzn.lambda[i]/(constants::c0/lBar*tBar));
        PIRobinBCValues[0] = -1.5*m_phtzn.lambda[i];
        PIRobinBCValues[1] = 1.5*m_phtzn.lambda[i];
        PIRobinBCValues[2] = 0;
        ParmParse PI("PI" + to_string(i+1));
        PI.getarr("bc_lo", PIBCLo, 0, SpaceDim);
        PI.getarr("bc_hi", PIBCHi, 0, SpaceDim);
      }
      else {
        for (int lev = startLev; lev<= finest_level; lev++)
          (s_PIAMRMG->levelOp(lev)).scale(*rhs[lev], m_phtzn.lambda[i]/m_phtzn.lambda[i-1]);
        ParmParse PI("PI" + to_string(i+1));
        PI.getarr("bc_lo", PIBCLo, 0, SpaceDim);
        PI.getarr("bc_hi", PIBCHi, 0, SpaceDim);
        PI.get("bc_value", PIDiriNeumBCValues);
      }
      
      int lbase = m_level;
      int lmax  = finest_level;
      s_PIAMRMG->solve(phi, rhs, lmax, lbase, false);
      
      for (int lev = m_level; lev <= finest_level; lev++)
        (*phi[lev]).copyTo(Interval(0,0), hierarchy[lev]->m_phtzn.Psi, Interval(i,i));
    }
    for (int lev = m_level; lev <= finest_level; lev++)
      hierarchy[lev]->m_phtzn.calcRate();
    
    //remember that startLev can be different from m_level
    for (int lev = startLev; lev<= finest_level; lev++) {
      delete rhs[lev];
      delete phi[lev];
      rhs[lev] = NULL;
      phi[lev] = NULL;
    }
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::tagCells " << m_level << endl;
  }
  
  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::tagCellsInit " << m_level << endl;
  }
  
  // Create tags based on undivided gradient of density
  IntVectSet localTags;
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
  {
    const AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
    
    PiecewiseLinearFillPatch pwl(levelDomain,
                                 amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
                                 1,
                                 amrGodCoarserPtr->m_problem_domain,
                                 amrGodCoarserPtr->m_ref_ratio,
                                 1);
    
    pwl.fillInterp(m_UNew,
                   amrGodCoarserPtr->m_UNew,
                   amrGodCoarserPtr->m_UNew,
                   1.0,
                   0,
                   0,
                   1);
  }
  m_UNew.exchange(Interval(0,1-1));
  
  // Compute undivided gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox gradFab(b,SpaceDim);
    const FArrayBox& UFab = m_UNew[dit()];
    
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const Box bCenter = b & grow(m_problem_domain,-BASISV(dir));
      const Box bLo     = b & adjCellLo(bCenter,dir);
      const int hasLo = ! bLo.isEmpty();
      const Box bHi     = b & adjCellHi(bCenter,dir);
      const int hasHi = ! bHi.isEmpty();
      FORT_GETGRADF(CHF_FRA1(gradFab,dir),
                    CHF_CONST_FRA1(UFab,0),
                    CHF_CONST_INT(dir),
                    CHF_BOX(bLo),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(bHi),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(bCenter));
    }
    
    FArrayBox gradMagFab(b,1);
    FORT_MAGNITUDEF(CHF_FRA1(gradMagFab,0),
                    CHF_CONST_FRA(gradFab),
                    CHF_BOX(b));
    
    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      
      if (gradMagFab(iv) >= m_refineThresh)
      {
        localTags |= iv;
      }
    }
  }
  
  localTags.grow(m_tagBufferSize);
  
  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;
  
  a_tags = localTags;
  
}

/*******/
void
AMRLevelAdvectDiffuse::
tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

/*******/
void
AMRLevelAdvectDiffuse::
regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3) {
    pout() << "AMRLevelAdvectDiffuse::regrid " << m_level << "; # of pts on this proc = " << m_grids.numPointsThisProc() << endl;
    printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "AMRLevelAdvectDiffuse::regrid-start");
  }
  
  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);
  
  // Save data for later
  LevelData<FArrayBox> UOld, ionOld;
  UOld.define(m_UNew);
  m_UNew.copyTo(m_UNew.interval(),UOld,UOld.interval());
  ionOld.define(m_ionNew);
  m_ionNew.copyTo(m_ionNew.interval(), ionOld, ionOld.interval());
  
  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,1,ivGhost);
  m_UOld.define(m_grids,1,ivGhost);
  m_dU.define(m_grids,1,ivGhost);
  m_ionNew.define(m_grids,m_gas.m_numOfIonSpe,ivGhost);
  m_ionOld.define(m_grids,m_gas.m_numOfIonSpe,ivGhost);
  
  LevelData<FArrayBox> neutOld;
  if (!m_gas.m_uniformity) {
    neutOld.define(m_neut);
    m_neut.copyTo(m_neut.interval(), neutOld, neutOld.interval());
    m_neut.define(m_grids,1,ivGhost);
  }
  
  IntVect ivGhost1 = m_numGhost*IntVect::Unit;
  LevelData<FArrayBox> phiOld;
  phiOld.define(m_phi);
  m_phi.copyTo(m_phi.interval(), phiOld, phiOld.interval());
    
  m_phi.define(m_grids, 1, ivGhost1);
  m_phiOld.define(m_grids, 1, ivGhost1);
  m_mu.define(m_grids, 1, ivGhost1);
  m_muOld.define(m_grids, 1, ivGhost1);
  m_advVel.define(m_grids,1,ivGhost);
  m_advVelOld.define(m_grids,1,ivGhost);
  
  m_field.define(m_grids, ivGhost1, ivGhost);
  m_fieldOld.define(m_grids, ivGhost1, ivGhost);
  
  LevelData<FluxBox> EEdgeOld, EEdgeCoarseOld;
  if (m_doImplicitPoisson) {
    s_aCoef[m_level] = (RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox> (m_grids, 1, IntVect::Zero)));
    s_bCoef[m_level] = (RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, IntVect::Zero)));
    s_aCoefComp[m_level] = (RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox> (m_grids, 1, IntVect::Zero)));
    s_bCoefComp[m_level] = (RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, IntVect::Unit)));
    EEdgeOld.define(m_EEdgeOld);
    EEdgeCoarseOld.define(m_EEdgeCoarseOld);
    m_EEdgeOld.copyTo(m_EEdgeOld.interval(), EEdgeOld, EEdgeOld.interval());
    m_EEdgeCoarseOld.copyTo(m_EEdgeCoarseOld.interval(), EEdgeCoarseOld, EEdgeCoarseOld.interval());
    m_EEdgeOld.define(m_grids,1,ivGhost);
    m_EEdgeCoarseOld.define(m_grids,1,ivGhost);
  }

  m_dUDiff.define(m_grids,1,ivGhost);
  m_flux.define(m_grids,1,ivGhost);
  m_muEdge.define(m_grids,1,ivGhost);
  
  photoionization PIOld;
  PIOld.define(m_phtzn.rate.disjointBoxLayout(), 3, ivGhost1);
  m_phtzn.copyTo(PIOld);

  m_phtzn.define(m_grids, 3, ivGhost1);
  m_phtznOld.define(m_grids, 3, ivGhost1);
  
  if (s_testing)
    m_testOutput.define(m_grids,1,ivGhost);
  
  // Set up data structures
  levelSetup();
  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    m_UNew[dit()].setVal(0.0);
    m_ionNew[dit()].setVal(0.0);
    if (!m_gas.m_uniformity)
      m_neut[dit()].setVal(0.0);
    m_phi[dit()].setVal(0.0);
    m_mu[dit()].setVal(0.0);
//    m_muOld[dit()].setVal(0.0);
    m_advVel[dit()].setVal(0.0);
    m_phtzn.Psi[dit()].setVal(0.0);
    m_phtzn.rate[dit()].setVal(0.0);
    m_flux[dit()].setVal(0.0);
    m_muEdge[dit()].setVal(0.0);
    if (m_doImplicitPoisson) {
      (*s_aCoef[m_level])[dit()].setVal(0.0);
      (*s_bCoef[m_level])[dit()].setVal(0.0);
      (*s_aCoefComp[m_level])[dit()].setVal(0.0);
      (*s_bCoefComp[m_level])[dit()].setVal(0.0);
      m_EEdgeOld[dit()].setVal(0.0);
      m_EEdgeCoarseOld[dit()].setVal(0.0);
    }
  }
  
  AMRLevelAdvectDiffuse* amrGodCoarserPtr;
  if (m_hasCoarser) amrGodCoarserPtr = getCoarserLevel();
  
  // Interpolate from coarser level
  if (m_hasCoarser) {
    //fills the valid region of current level by piecewise linear interpolation from data on a coarser level of refinement,
    m_fineInterp.interpToFine(m_UNew, amrGodCoarserPtr->m_UNew);
    m_fineInterpM.interpToFine(m_ionNew, amrGodCoarserPtr->m_ionNew);
    if (!m_gas.m_uniformity)
      m_fineInterp.interpToFine(m_neut, amrGodCoarserPtr->m_neut);
    m_fineInterp.interpToFine(m_phi, amrGodCoarserPtr->m_phi);
    m_fineInterpPhtzn.interpToFine(m_phtzn.Psi, amrGodCoarserPtr->m_phtzn.Psi);
    // revisit: need interpolation in time?
    if (m_doImplicitPoisson) {
      m_fineInterpFace.interpToFine(m_EEdgeOld, amrGodCoarserPtr->m_EEdgeOld);
      m_fineInterpFace.interpToFine(m_EEdgeCoarseOld, amrGodCoarserPtr->m_EEdgeOld);
    }
  }
  
  if (m_hasCoarser) {
    Real eps = 1.0e-10;
    // need this otherwise m_xxxOld will have invalid values when being assigned by copying at the beginning of advance, because fineInterp doesn't assign values to ghost cells. Otherwise, we will have problems if the ghost cells of m_xxxOld are used somewhere in the code. For example, in fillGhost with time interpolation.
    fillGhost(m_pwl, m_UNew, amrGodCoarserPtr->m_UNew, 0, amrGodCoarserPtr->m_UNew, 0, eps, 0);
    fillGhost(m_pwlM, m_ionNew, amrGodCoarserPtr->m_ionNew, 0, amrGodCoarserPtr->m_ionNew, 0, eps, 0);
    if (!m_gas.m_uniformity)
      fillGhost(m_pwl, m_neut, amrGodCoarserPtr->m_neut, 0, amrGodCoarserPtr->m_neut, 0, eps, 0);
    fillGhost(m_pwl, m_phi, amrGodCoarserPtr->m_phi, 0, amrGodCoarserPtr->m_phi, 0, eps, 0);
  }
  // Copy from old state
  UOld.copyTo(UOld.interval(), m_UNew, m_UNew.interval());
  ionOld.copyTo(ionOld.interval(), m_ionNew, m_ionNew.interval());
  phiOld.copyTo(phiOld.interval(), m_phi, m_phi.interval());
  PIOld.copyTo(m_phtzn);
  
  if (!m_gas.m_uniformity)
    neutOld.copyTo(neutOld.interval(), m_neut, m_neut.interval());
  
  if (m_doImplicitPoisson) {
    EEdgeOld.copyTo(EEdgeOld.interval(), m_EEdgeOld, m_EEdgeOld.interval());
    EEdgeCoarseOld.copyTo(EEdgeCoarseOld.interval(), m_EEdgeCoarseOld, m_EEdgeCoarseOld.interval());
  }
//  m_phtzn.calcRate();
  
  computeEField(false);
  fillMobility(false);
  fillAdvectionVelocity(false);
  if (m_phtzn.runSolve) {
    m_phtzn.calcRate();
  }
  if (s_verbosity >= 3) {
    pout() << "AMRLevelAdvectDiffuse::regrid-end " << m_level << " pts of this proc = " << m_grids.numPointsThisProc() << endl;
    printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "AMRLevelAdvectDiffuse::regrid-end");
  }
}


/*******/
void
AMRLevelAdvectDiffuse::
initialGrid(const Vector<Box>& a_newGrids)
{
  
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::initialGrid " << m_level << endl;
  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);
  
  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,1,ivGhost);
  m_UOld.define(m_grids,1,ivGhost);
  m_dU.define(m_grids,1,ivGhost);
  m_ionNew.define(m_grids,m_gas.m_numOfIonSpe,ivGhost);
  m_ionOld.define(m_grids,m_gas.m_numOfIonSpe,ivGhost);
  if (!m_gas.m_uniformity)
    m_neut.define(m_grids,1,ivGhost);
  
  IntVect ivGhost1 = m_numGhost*IntVect::Unit;
  m_phi.define(m_grids, 1, ivGhost1);
  m_phiOld.define(m_grids, 1, ivGhost1);
  m_mu.define(m_grids, 1, ivGhost1);
  m_muOld.define(m_grids, 1, ivGhost1);
  m_advVel.define(m_grids,1,ivGhost);
  m_advVelOld.define(m_grids,1,ivGhost);
  
  m_dUDiff.define(m_grids,1,ivGhost);
  m_flux.define(m_grids,1,ivGhost);
  m_muEdge.define(m_grids,1,ivGhost);
  
  m_field.define(m_grids, ivGhost1, ivGhost);
  m_fieldOld.define(m_field);
  
  if (m_doImplicitPoisson) {
    s_aCoef[m_level] = (RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox> (m_grids, 1, IntVect::Zero)));
    s_bCoef[m_level] = (RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, IntVect::Zero)));
    s_aCoefComp[m_level] = (RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox> (m_grids, 1, IntVect::Zero)));
    s_bCoefComp[m_level] = (RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, IntVect::Unit)));
    m_EEdgeOld.define(m_grids,1,ivGhost);
    m_EEdgeCoarseOld.define(m_grids,1,ivGhost);
  }
  
  m_phtzn.define(m_grids, 3, ivGhost1);
  m_phtznOld.define(m_grids, 3, ivGhost1);
  
  if (s_testing)
    m_testOutput.define(m_grids,1,ivGhost);
  
  // Set up data structures
  levelSetup();
}

/*******/
void
AMRLevelAdvectDiffuse::
initialData()
{
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::initialData " << m_level << endl;

  PhysIBC* physIBCPtr = m_advPhys->getPhysIBC();
  physIBCPtr->initialize(m_UNew);
    
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    const Box& b = m_UNew[dit()].box();
    m_UOld[dit()].copy(m_UNew[dit()], b);
    // positive ions
    m_ionNew[dit()].copy(m_UNew[dit()], b, 0, b, 0, 1);
    // negative ions
    m_ionNew[dit()].setVal(0.0, 1);
    m_ionOld[dit()].copy(m_ionNew[dit()], b);
    m_phi[dit()].setVal(0.0, m_phi[dit()].box(), 0);
    m_mu[dit()].setVal(0.0, m_mu[dit()].box(), 0);
    m_muOld[dit()].setVal(0.0, m_muOld[dit()].box(), 0);
    m_field.m_E[dit()].setVal(0.0, m_field.m_E[dit()].box(), 0, SpaceDim);
    m_field.m_Emag[dit()].setVal(0.0, m_field.m_Emag[dit()].box(), 0);
    m_field.m_EEdge[dit()].setVal(0.0, m_field.m_EEdge[dit()].box());
    m_advVel[dit()].setVal(0.0, m_advVel[dit()].box());
    m_advVelOld[dit()].setVal(0.0, m_advVelOld[dit()].box());
    m_phtzn.Psi[dit()].setVal(0.0, m_phtzn.Psi[dit()].box(), 0, 3);
    m_phtzn.rate[dit()].setVal(0.0, m_phtzn.rate[dit()].box(), 0);
    m_flux[dit()].setVal(0.0);
    m_muEdge[dit()].setVal(0.0);
    
    if (m_doImplicitPoisson) {
      (*s_aCoef[m_level])[dit()].setVal(0.0);
      (*s_bCoef[m_level])[dit()].setVal(0.0);
      (*s_aCoefComp[m_level])[dit()].setVal(0.0);
      (*s_bCoefComp[m_level])[dit()].setVal(0.0);
      m_EEdgeOld[dit()].setVal(0.0);
      m_EEdgeCoarseOld[dit()].setVal(0.0);
    }
    
    if (s_testing)
      m_testOutput[dit()].setVal(0.0);
  }
  
  ParmParse pp("bgdChargeCloud");
  if (pp.contains("number")) {
    int num;
    pp.get("number", num);
        
    if (!pp.contains("center")) {
      if (m_level == 0) {
        Vector<Real> center(SpaceDim*num, 0.0);
        Vector<Real> distCenter(SpaceDim, 0.0);
        Vector<Real> distBoxLength(SpaceDim, 0.0);
        
        pp.getarr("distCenter", distCenter, 0, SpaceDim);
        pp.getarr("distBoxLength", distBoxLength, 0, SpaceDim);
        
        const std::string varName("bgdChargeCloud.center");
        std::string valStr;
        if (procID() == uniqueProc(SerialTask::compute))
          for (int i = 0; i < num; i++)
            for (int dir = 0; dir < SpaceDim; dir++) {
              center[dir+i*SpaceDim] = ((1.0*rand()/RAND_MAX) - 0.5) * distBoxLength[dir] + distCenter[dir];
              valStr += to_string(center[dir+i*SpaceDim]);
              valStr += " ";
            }
        broadcast(valStr,uniqueProc(SerialTask::compute));
        pp.setStr(varName, valStr);
      }
    }
    
    Real r0, mag0;
    Vector<Real> center(SpaceDim*num, 0.0);
    pp.get("radius", r0);
    r0 /= lBar;
    pp.get("mag", mag0);
    mag0 /= nBar;
    pp.getarr("center", center, 0, SpaceDim*num);
    for (int i = 0; i < SpaceDim*num; i++)
     center[i] /= lBar;
    
    Vector<Real> radius(SpaceDim*num, r0), mag(SpaceDim*num, mag0);
//    pp.getarr("radius", radius, 0, SpaceDim*num);
//    pp.getarr("mag",   mag, 0, num);
    
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
      for (BoxIterator bit(m_grids.get(dit)); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        for (int i = 0; i < num; i++) {
          double val = fabs(mag[i]);
          for (int dir = 0; dir < SpaceDim; dir++) {
            double tmp = 0;
            tmp = (iv[dir]+0.5)*m_dx - center[dir+i*SpaceDim];
            tmp = tmp*tmp/radius[dir+i*SpaceDim]/radius[dir+i*SpaceDim];
            val *= exp(-tmp);
          }
          if (mag[i] > 0)
            m_ionNew[dit()](iv, 0) += val;
          else
            m_ionNew[dit()](iv, 1) += val;
        }
      }
      m_ionOld[dit()].copy(m_ionNew[dit()], m_ionOld[dit()].box());
    }
  }
  
  /* for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    const Box& b = m_UNew[dit()].box();
    
    if (!m_gas.m_uniformity)
      m_neut[dit()].setVal(m_gas.m_N, b, 0);
  } */
  if (!m_gas.m_uniformity) {
    Real levDx = m_dx;
    RealVect ccOffset = 0.5*levDx*RealVect::Unit;
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
      for (BoxIterator bit(m_grids.get(dit)); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        
        RealVect point(iv);
        point *= levDx;
        point += ccOffset;
/*          point *= normalization::lBar;
        m_neut[dit()](iv, 0) = neutDensityFile->value(point) / normalization::nBar; */
        vector< double> pvec(SpaceDim);
        for(int i = 0; i != pvec.size(); i++)
          pvec[i] = point[i];
        m_neut[dit()](iv, 0) = m_gas.getBackgroundDensity(pvec);
      }
    }
  }
  
  pp = ParmParse("gas");
  if (pp.contains("inhom")) {
    bool flag;
    pp.get("inhom", flag);
    
    if (flag) {
      if (pp.contains("inhomRandom")) {
        bool flagForRandom;
        pp.get("inhomRandom", flagForRandom);
        if (flagForRandom) {
          double inhomDensity;
          pp.get("inhomDensity", inhomDensity);
          inhomDensity *= lBar*lBar*lBar;
          
          IntVect sizeVec = m_problem_domain.domainBox().size();
          
          for (int idir = 0; idir < SpaceDim; idir++)
            inhomDensity *= m_dx * sizeVec[idir];
          int num = inhomDensity;
          
          Real r0, mag0;
          pp.get("inhomRadius", r0);
          pp.get("inhomMag", mag0);
          r0 /= lBar;
          Vector<Real> radius(SpaceDim*num, r0), mag(SpaceDim*num, mag0);
          
          Vector<Real> center(SpaceDim*num, 0.0);          
          if (!pp.contains("inhomCenter")) {
            Vector<Real> center(SpaceDim*num, 0.0);
            const std::string varName("gas.inhomCenter");
            std::string valStr;
            if (procID() == uniqueProc(SerialTask::compute)) {
              for (int i = 0; i < num; i++) {
                for (int dir = 0; dir < SpaceDim; dir++) {
                  center[dir+i*SpaceDim] = (1.0*rand()/RAND_MAX) * sizeVec[dir] * m_dx;
                  valStr += to_string(center[dir+i*SpaceDim]);
                  valStr += " ";
                }
              }
            }
            broadcast(valStr, uniqueProc(SerialTask::compute));
            broadcast(center,uniqueProc(SerialTask::compute));
            pp.setStr(varName, valStr);
          } else {
            pp.getarr("inhomCenter", center, 0, SpaceDim*num);
          }
          
          for (int i = 0; i < num; i++) {
            vector< double> pvec(SpaceDim);
            for(int dir = 0; dir != pvec.size(); dir++) {
              pvec[dir] = center[dir+i*SpaceDim];
            }
            double val0 = mag[i] * m_gas.getBackgroundDensity(pvec);
            
            for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
              for (BoxIterator bit(m_grids.get(dit)); bit.ok(); ++bit) {
                const IntVect& iv = bit();
                double val = val0;
                for (int dir = 0; dir < SpaceDim; dir++) {
                  double tmp = 0;
                  tmp = (iv[dir]+0.5)*m_dx - center[dir+i*SpaceDim];
                  tmp = tmp*tmp/radius[dir+i*SpaceDim]/radius[dir+i*SpaceDim];
                  val *= exp(-tmp);
                }
                m_neut[dit()](iv, 0) += val;
              }
            }
          }
        }
      }
    }
  }
  
  m_phi.copyTo(m_phiOld);
  m_phtzn.copyTo(m_phtznOld);
  
}

/*******/
void
AMRLevelAdvectDiffuse::
postInitialize()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::postInitialize " << m_level << endl;
  }
  
  // poissonSolve is done at m_level = 0 at this step because time at all levels is 0;
  poissonSolveComposite();
  
  if (s_verbosity >= 3) {
    printDiagnosticInfo (m_level, m_dx, m_grids, m_UNew, "U", "postInitialize");
//    outputDataForCheck (m_level, m_grids, m_UNew);
//    printDiagnosticInfo (m_level, m_dx, m_grids, m_ionNew, "ion", "postInitialize");
//    printDiagnosticInfo (m_level, m_dx, m_grids, m_advVel, "advVel", "postInitialize");
//    outputDataForCheck (m_level, m_grids, m_phi);
  }
  
  if (m_level == 0) {
    
    Vector<AMRLevelAdvectDiffuse*>         hierarchy;
    Vector<int>                            refRat;
    Vector<DisjointBoxLayout>              grids;
    Real                                   lev0Dx;
    ProblemDomain                          lev0Domain;
    
    getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
    int finest_level = hierarchy.size()-1;
    for (int lev= m_level; lev <= finest_level; lev++) {
      hierarchy[lev]->fillMobility(false);
      hierarchy[lev]->fillAdvectionVelocity(false);
      
      if (m_doImplicitPoisson) {
        hierarchy[lev]->m_field.m_EEdge.copyTo(hierarchy[lev]->m_EEdgeOld);
        if (hierarchy[lev]->m_hasFiner)
          hierarchy[lev+1]->m_field.m_EEdge.copyTo(hierarchy[lev+1]->m_EEdgeCoarseOld);
      }
      
      if (s_testing)
        hierarchy[lev]->testing();
    }
  }
  if (m_phtzn.runSolve) {
    photoionizationSolve();
  }
  
  outputStepStats(AMPLIFIOut);
}

#ifdef CH_USE_HDF5

/*******/
void
AMRLevelAdvectDiffuse::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::writeCheckpointHeader" << endl;
  }
  
  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = 1;
  
  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < 1; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }
  
  // Write the header
  header.writeToFile(a_handle);
  
  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::writeCheckpointLevel" << endl;
  }
  
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;
  
  a_handle.setGroup(label);
  
  // Setup the level header information
  HDF5HeaderData header;
  
  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = m_dx;
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_real["nu"]              = m_nu;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();
  
  
  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
         header.m_int ["is_periodic_0"] = 1;
         else
         header.m_int ["is_periodic_0"] = 0; ,
         
         if (m_problem_domain.isPeriodic(1))
         header.m_int ["is_periodic_1"] = 1;
         else
         header.m_int ["is_periodic_1"] = 0; ,
         
         if (m_problem_domain.isPeriodic(2))
         header.m_int ["is_periodic_2"] = 1;
         else
         header.m_int ["is_periodic_2"] = 0; );
  
  // Write the header for this level
  header.writeToFile(a_handle);
  
  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
  
  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data");
}

/*******/
void
AMRLevelAdvectDiffuse::
readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::readCheckpointHeader" << endl;
  }
  
  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  
  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }
  
  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: checkpoint file does not have num_components");
  }
  
  int numStates = header.m_int["num_components"];
  if (numStates != 1)
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }
  
  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < 1; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: checkpoint file does not have enough component names");
    }
    
    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::readCheckpointLevel" << endl;
  }
  
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;
  
  // Read the header for this level
  a_handle.setGroup(label);
  
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  
  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }
  
  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];
  
  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }
  
  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];
  
  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }
  
  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];
  
  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }
  
  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];
  
  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }
  
  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];
  
  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }
  
  // Get nu
  if (header.m_real.find("nu") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain nu");
  }
  m_nu = header.m_real["nu"];
  
  if (s_verbosity >= 2)
  {
    pout() << "read nu = " << m_nu << endl;
  }
  
  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain prob_domain");
  }
  
  Box domainBox = header.m_box["prob_domain"];
  
  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
         isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
         isPeriodic[0] = false; ,
         
         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
         isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
         isPeriodic[1] = false; ,
         
         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
         isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
         isPeriodic[2] = false;);
  
  m_problem_domain = ProblemDomain(domainBox,isPeriodic);
  
  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);
  
  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain a Vector<Box>");
  }
  
  // Create level domain
  Vector<int> procs;
  LoadBalance(procs, grids);
  m_grids = DisjointBoxLayout(grids,procs, m_problem_domain);
  
  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;
  
  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = constGrids[lit()];
    m_level_grids.push_back(b);
  }
  
  if (s_verbosity >= 4)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = m_grids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }
  
  // Reshape state with new grids
  m_UNew.define(m_grids,1);
  m_dU.define(  m_grids,1);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_UNew,
                                         "data",
                                         m_grids);
  
  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain state data");
  }
  m_UOld.define(m_grids,1);
  m_advVel.define(m_grids,1, m_numGhost*IntVect::Unit);
  fillMobility(false);
  fillAdvectionVelocity(false);
  
  // Set up data structures
  levelSetup();
}

/*******/
void
AMRLevelAdvectDiffuse::
writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::writePlotHeader" << endl;
  }

  // Setup the number of components -- include space for error
  HDF5HeaderData header;

  header.m_real["time"] = m_time * normalization::tBar;
  
  // Setup the component names
  char compStr[30];
  int istart = 0;

  for (int comp = 0; comp < m_UNew.nComp(); ++comp)
  {
    sprintf(compStr,"component_%d",istart+comp);
    header.m_string[compStr] = m_stateNames[comp];
  }
  istart += m_UNew.nComp();

  for (int comp = 0; comp < m_ionNew.nComp(); ++comp) {
    sprintf(compStr,"component_%d",istart+comp);
    header.m_string[compStr] = "ion"+to_string(comp);
  }
  istart += m_ionNew.nComp();
  
  if (!m_gas.m_uniformity) {
    sprintf(compStr,"component_%d",istart);
    header.m_string[compStr] = "neut";
    istart += m_neut.nComp();
  }
  
  sprintf(compStr,"component_%d",istart);
  header.m_string[compStr] = "E";
  istart += m_field.m_Emag.nComp();
  
//  sprintf(compStr,"component_%d",istart);
//  header.m_string[compStr] = "EOld";
//  istart += m_fieldOld.m_Emag.nComp();
  
//  sprintf(compStr,"component_%d",istart);
//  header.m_string[compStr] = "srs";
//  istart += m_srs.nComp();
  
  sprintf(compStr,"component_%d",istart);
  header.m_string[compStr] = "PI";
  istart += m_phtzn.rate.nComp();
  
//  sprintf(compStr,"component_%d",istart);
//  for (int comp = 0; comp < m_phtzn.Psi.nComp(); ++comp) {
//    sprintf(compStr,"component_%d",istart+comp);
//    header.m_string[compStr] = "PI"+to_string(comp+1);
//  }
//  istart += m_phtzn.Psi.nComp();

  sprintf(compStr,"component_%d",istart);
  header.m_string[compStr] = "phi";
  istart += m_phi.nComp();

//  sprintf(compStr,"component_%d",istart);
//  header.m_string[compStr] = "phiOld";
//  istart += m_phiOld.nComp();
//
//  for (int comp = 0; comp < m_field.m_E.nComp(); ++comp) {
//    sprintf(compStr,"component_%d",istart+comp);
//    header.m_string[compStr] = "E"+to_string(comp+1);
//  }
//  istart += m_field.m_E.nComp();
//
//  for (int comp = 0; comp < m_fieldOld.m_E.nComp(); ++comp) {
//    sprintf(compStr,"component_%d",istart+comp);
//    header.m_string[compStr] = "EOld"+to_string(comp+1);
//  }
//  istart += m_fieldOld.m_E.nComp();

//  if(m_doImplicitPoisson) {
//    for (int dir=0; dir<SpaceDim; dir++) {
//      sprintf(compStr,"component_%d",istart+dir);
//      header.m_string[compStr] = "EEdgeToCell"+to_string(dir);
//    }
//    istart += m_field.m_EEdge.nComp()*SpaceDim;
//
//    for (int dir=0; dir<SpaceDim; dir++) {
//      sprintf(compStr,"component_%d",istart+dir);
//      header.m_string[compStr] = "EEdgeToCellOld"+to_string(dir);
//    }
//    istart += m_fieldOld.m_EEdge.nComp()*SpaceDim;
//  }
  
//  sprintf(compStr,"component_%d",istart);
//  header.m_string[compStr] = "mu";
//  istart += m_mu.nComp();
  
//  for (int dir=0; dir<SpaceDim; dir++) {
//    sprintf(compStr,"component_%d",istart+dir);
//    header.m_string[compStr] = "VeEdgeToCell"+to_string(dir);
//  }
//  istart += m_advVel.nComp()*SpaceDim;
  
//  for (int dir=0; dir<SpaceDim; dir++) {
//    sprintf(compStr,"component_%d",istart+dir);
//    header.m_string[compStr] = "flux"+to_string(dir);
//  }
//  istart += m_flux.nComp()*SpaceDim;
  
//  if(m_doImplicitPoisson) {
//    for (int dir=0; dir<SpaceDim; dir++) {
//      sprintf(compStr,"component_%d",istart+dir);
//      header.m_string[compStr] = "bCoeff"+to_string(dir);
//    }
//    istart += (*s_bCoef[m_level]).nComp()*SpaceDim;
//  }
//
//  if(m_doImplicitPoisson) {
//    for (int dir=0; dir<SpaceDim; dir++) {
//      sprintf(compStr,"component_%d",istart+dir);
//      header.m_string[compStr] = "bCoeffComp"+to_string(dir);
//    }
//    istart += (*s_bCoefComp[m_level]).nComp()*SpaceDim;
//  }
  
  if (s_testing) {
    sprintf(compStr,"component_%d",istart);
    header.m_string[compStr] = "test";
    istart += m_testOutput.nComp();
  }
    
  header.m_int["num_components"] = istart;

  // Write the header
  header.writeToFile(a_handle);

  //  char levelStr[20];
  //  sprintf(levelStr,"%d",m_level);
  //  const std::string label = std::string("level_") + levelStr;
  //  a_handle.setGroup(label+"flux");
  //  header.m_int["num_components"] = 1;
  //  sprintf(compStr,"component_%d", 0);
  //  header.m_string[compStr] = "EEdge";
  //  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::writePlotLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx * normalization::lBar;
  header.m_real["dt"]          = m_dt * normalization::tBar;
  header.m_real["time"]        = m_time * normalization::tBar;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  int numPlotVar;
  
//  int numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + m_field.m_E.nComp();

//  int numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + m_field.m_E.nComp() + m_advVel.nComp()*SpaceDim;
  
  
  
//  int numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phi.nComp() + m_field.m_E.nComp() + m_field.m_EEdge.nComp()*SpaceDim;
  if (m_doImplicitPoisson)
    if (!m_gas.m_uniformity)
      numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_neut.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + s_testing*m_testOutput.nComp();
    else
      numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + s_testing*m_testOutput.nComp();
//    numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + m_phiOld.nComp() + m_field.m_E.nComp() + m_fieldOld.m_E.nComp() + (*s_bCoef[m_level]).nComp()*SpaceDim + (*s_bCoefComp[m_level]).nComp()*SpaceDim;
  else
//    numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + m_flux.nComp()*SpaceDim;
    if (!m_gas.m_uniformity)
      numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_neut.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + s_testing*m_testOutput.nComp();
    else
      numPlotVar = m_UNew.nComp() + m_ionNew.nComp() + m_field.m_Emag.nComp() + m_phtzn.rate.nComp() + m_phi.nComp() + s_testing*m_testOutput.nComp();
  
//  numPlotVar += m_mu.nComp();
  
  LevelData<FArrayBox> plotData(m_UNew.disjointBoxLayout(), numPlotVar);

  // first copy data to plot data holder
  int istart = 0, nComp;
  nComp = m_UNew.nComp();
  Interval interv(istart, istart+nComp-1);
  m_UNew.copyTo(m_UNew.interval(), plotData, interv);
  unnormalize(plotData, istart, nComp, normalization::nBar);
  istart += nComp;

  nComp = m_ionNew.nComp();
  interv.define(istart, istart+nComp-1);
  m_ionNew.copyTo(m_ionNew.interval(), plotData, interv);
  unnormalize(plotData, istart, nComp, normalization::nBar);
  istart += nComp;

  if (!m_gas.m_uniformity) {
    nComp = m_neut.nComp();
    interv.define(istart, istart+nComp-1);
    m_neut.copyTo(m_neut.interval(), plotData, interv);
    unnormalize(plotData, istart, nComp, normalization::nBar);
    istart += nComp;
  }
  
  nComp = m_field.m_Emag.nComp();
  interv.define(istart, istart+nComp-1);
  m_field.m_Emag.copyTo(m_field.m_Emag.interval(), plotData, interv);
  unnormalize(plotData, istart, nComp, normalization::EBar);
  istart += nComp;
  
//  nComp = m_fieldOld.m_Emag.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_fieldOld.m_Emag.copyTo(m_fieldOld.m_Emag.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::EBar);
//  istart += nComp;

//  nComp = m_srs.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_srs.copyTo(m_srs.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, 1/(pow(normalization::lBar,3)*normalization::tBar));
//  istart += nComp;
  
  nComp = m_phtzn.rate.nComp();
  interv.define(istart, istart+nComp-1);
  m_phtzn.rate.copyTo(m_phtzn.rate.interval(), plotData, interv);
  unnormalize(plotData, istart, nComp, 1/(pow(normalization::lBar,3)*normalization::tBar));
  istart += nComp;
  
//  nComp = m_phtzn.Psi.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_phtzn.Psi.copyTo(m_phtzn.Psi.interval(), plotData, interv);
//  for (int comp = 0; comp <m_phtzn.Psi.nComp(); comp++)
//    unnormalize(plotData, istart+comp, 1, m_phtzn.A[comp]/pow(lBar,4)*constants::c0);
//  istart += nComp;
  
  nComp = m_phi.nComp();
  interv.define(istart, istart+nComp-1);
  m_phi.copyTo(m_phi.interval(), plotData, interv);
  unnormalize(plotData, istart, nComp, normalization::phiBar);
  istart += nComp;
  
//  nComp = m_phiOld.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_phiOld.copyTo(m_phiOld.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::phiBar);
//  istart += nComp;
//
//  nComp = m_field.m_E.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_field.m_E.copyTo(m_field.m_E.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::EBar);
//  istart += nComp;
//
//  nComp = m_fieldOld.m_E.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_fieldOld.m_E.copyTo(m_fieldOld.m_E.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::EBar);
//  istart += nComp;

//  if (m_doImplicitPoisson) {
//    LevelData<FArrayBox> tmp1(m_UNew.getBoxes(), SpaceDim);
//    EdgeToCell(m_field.m_EEdge, tmp1);
//    nComp = tmp1.nComp();
//    interv.define(istart, istart+nComp-1);
//    tmp1.copyTo(tmp1.interval(), plotData, interv);
//    unnormalize(plotData, istart, nComp, normalization::EBar);
//    istart += nComp;
//
//    EdgeToCell(m_fieldOld.m_EEdge, tmp1);
//    nComp = tmp1.nComp();
//    interv.define(istart, istart+nComp-1);
//    tmp1.copyTo(tmp1.interval(), plotData, interv);
//    unnormalize(plotData, istart, nComp, normalization::EBar);
//    istart += nComp;
//  }
  
//  nComp = m_mu.nComp();
//  interv.define(istart, istart+nComp-1);
//  m_mu.copyTo(m_mu.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::muBar);
//  istart += nComp;
  
//  nComp = tmp.nComp();
//  interv.define(istart, istart+nComp-1);
//  tmp.copyTo(tmp.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::EBar);
//  istart += nComp;
  
//  LevelData<FArrayBox> tmp(m_UNew.getBoxes(), SpaceDim);
//  EdgeToCell(m_advVel, tmp);
//  nComp = tmp.nComp();
//  interv.define(istart, istart+nComp-1);
//  tmp.copyTo(tmp.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::lBar/normalization::tBar);
//  istart += nComp;
  
//  LevelData<FArrayBox> tmp(m_UNew.getBoxes(), SpaceDim);
//  EdgeToCell(m_flux, tmp);
//  nComp = tmp.nComp();
//  interv.define(istart, istart+nComp-1);
//  tmp.copyTo(tmp.interval(), plotData, interv);
//  unnormalize(plotData, istart, nComp, normalization::nBar*lBar/tBar);
//  istart += nComp;

//  LevelData<FArrayBox> tmp(m_UNew.getBoxes(), SpaceDim);
//  if (m_doImplicitPoisson) {
//    EdgeToCell((*s_bCoef[m_level]), tmp);
//    nComp = tmp.nComp();
//    interv.define(istart, istart+nComp-1);
//    tmp.copyTo(tmp.interval(), plotData, interv);
//  //  unnormalize(plotData, istart, nComp, normalization::nBar*lBar/tBar);
//    unnormalize(plotData, istart, nComp, 1.0);
//    istart += nComp;
//  }
//
//  if (m_doImplicitPoisson) {
//    EdgeToCell((*s_bCoefComp[m_level]), tmp);
//    nComp = tmp.nComp();
//    interv.define(istart, istart+nComp-1);
//    tmp.copyTo(tmp.interval(), plotData, interv);
//  //  unnormalize(plotData, istart, nComp, normalization::nBar*lBar/tBar);
//    unnormalize(plotData, istart, nComp, 1.0);
//    istart += nComp;
//  }
  
  if (s_testing) {
    nComp = m_testOutput.nComp();
    interv.define(istart, istart+nComp-1);
    m_testOutput.copyTo(m_testOutput.interval(), plotData, interv);
    unnormalize(plotData, istart, nComp, 1.0);
    istart += nComp;
  }
  write(a_handle,m_UNew.disjointBoxLayout());
  write(a_handle,plotData,"data");

  //  a_handle.setGroup(label+"flux");
  //  header.writeToFile(a_handle);
  //  LevelData<FluxBox> plotDataFlux(m_field.m_EEdge.getBoxes(), 1);
  //  m_field.m_EEdge.copyTo(m_field.m_EEdge.interval(), plotDataFlux, m_field.m_EEdge.interval());
  //  write(a_handle,m_field.m_EEdge.boxLayout());
  //  write(a_handle,plotDataFlux,"flux");
  
//  a_handle.setGroup(std::string("flux_level_") + levelStr);
//  header.writeToFile(a_handle);
//  LevelData<FluxBox> plotDataFlux(m_advVel.getBoxes(), 1);
//  m_advVel.copyTo(plotDataFlux);
//  write(a_handle,m_advVel.boxLayout());
//  write(a_handle,plotDataFlux,"flux");
}

#endif

/*******/
Real
AMRLevelAdvectDiffuse::
computeDt()
{
  
  Real newDtC, eps = 1e-20;
  
  if(m_levelGodunov.getMaxWaveSpeed(m_UNew, m_advVel) > eps)
    newDtC = m_cfl * m_dx / m_levelGodunov.getMaxWaveSpeed(m_UNew, m_advVel);
  else if (m_nu > eps)
    newDtC = m_cfl*m_dx*m_dx/m_nu;
  else
    MayDay::Error("Insufficient condition to determine time step");
    
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::computeDt " << m_level << " max speed " << m_levelGodunov.getMaxWaveSpeed(m_UNew, m_advVel) << " dtC " << newDtC*tBar << endl;
  
  Real newDT;
  newDT = min(newDtC, computeDtI());
  newDT = min(newDT, computeDtM());
  
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::computeDt " << m_level << " new dt " << newDT*tBar << endl;
  
  return newDT;
}

///*******/
//Real
//AMRLevelAdvectDiffuse::computeDtI()
//{
//  Real Emax = 0;
//
//  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
//    Emax = max(m_field.m_Emag[dit()].max(m_grids[dit()], 0), Emax);
//
//  // Gather maximum wave speeds and broadcast the maximum over these
//  Vector<Real> allEmax;
//  gather(allEmax,Emax,uniqueProc(SerialTask::compute));
//
//  Real eps = 1e-20;
//  Real dtI;
//  if (procID() == uniqueProc(SerialTask::compute)) {
//    Real ai, ve;
//    Emax = allEmax[0];
//    for (int i = 1; i < allEmax.size (); ++i)
//      Emax = max(Emax,allEmax[i]);
//
//    Real n = m_gas.m_N;
//
//    ai = getProcessCoeff(Emax/n, m_gas, "ionization")*n;
//    ve = getProcessCoeff(Emax/n, m_gas, "mobility")/n*Emax;
//    dtI = numerical::Ai * 1/max(ai*ve, eps);
//  }
//
//  broadcast(dtI,uniqueProc(SerialTask::compute));
//  broadcast(Emax,uniqueProc(SerialTask::compute));
//  if (s_verbosity >= 3)
//    pout() << "AMRLevelAdvectDiffuse::computeDtI " << m_level << " Emax " << Emax*EBar << " dtI " << dtI*tBar << endl;
//  return dtI;
//}

/*******/
Real
AMRLevelAdvectDiffuse::computeDtI()
{
  if (m_gas.m_uniformity) {
    Real Emax = 0;
       
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      Emax = max(m_field.m_Emag[dit()].max(m_grids[dit()], 0), Emax);

    // Gather maximum Emax and broadcast the maximum over these
    Vector<Real> allEmax;
    gather(allEmax,Emax,uniqueProc(SerialTask::compute));
    
    Real eps = 1e-20;
    Real dtI;
    if (procID() == uniqueProc(SerialTask::compute)) {
      Real ai, ve;
      Emax = allEmax[0];
      for (int i = 1; i < allEmax.size(); ++i)
        Emax = max(Emax,allEmax[i]);
      
      Real n = m_gas.m_N;
      
      ai = m_gas.EDpdentProcs["ionization"].value(Emax/n)*n;
      ve = m_gas.EDpdentProcs["mobility"].value(Emax/n)/n*Emax;
      dtI = numerical::Ai * 1/max(ai*ve, eps);
    }
    broadcast(dtI,uniqueProc(SerialTask::compute));
    broadcast(Emax,uniqueProc(SerialTask::compute));
    if (s_verbosity >= 3)
      pout() << "AMRLevelAdvectDiffuse::computeDtI " << m_level << " Emax " << Emax*EBar << " dtI " << dtI*tBar << endl;
    return dtI;
  }
  
  else {
    Real rateImax = 0;
    
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
      FArrayBox rateI;
      rateI.define(m_field.m_Emag[dit()].box(), 1);
      getRate(rateI, m_field.m_Emag[dit()], m_neut[dit()], m_mu[dit()], "ionization");
      rateImax = max(rateI.max(m_grids[dit()], 0), rateImax);
    }
    // Gather maximum rate and broadcast the maximum over these
    Vector<Real> allRateImax;
    gather(allRateImax, rateImax, uniqueProc(SerialTask::compute));
    
    Real eps = 1e-20;
    Real dtI;
    if (procID() == uniqueProc(SerialTask::compute)) {
      rateImax = allRateImax[0];
      for (int i = 1; i < allRateImax.size(); ++i)
        rateImax = max(rateImax, allRateImax[i]);
      dtI = numerical::Ai * 1/max(rateImax, eps);
    }
    broadcast(dtI,uniqueProc(SerialTask::compute));
    if (s_verbosity >= 3)
      pout() << "AMRLevelAdvectDiffuse::computeDtI " << m_level << " dtI " << dtI*tBar << endl;
    return dtI;
  }
}

/*******/
Real
AMRLevelAdvectDiffuse::computeDtM()
{
  DataIterator dit = m_grids.dataIterator();

  Real dtM = numeric_limits<double>::min();
  
  LevelData<FArrayBox> muU;
  muU.define(m_grids, 1, m_numGhost*IntVect::Unit);
  m_UNew.copyTo(m_UNew.interval(), muU, muU.interval());
  // Loop over all grids
  for (dit.begin(); dit.ok(); ++dit) {
    const Box& curBox = m_grids[dit()];

    FArrayBox& muFab = m_mu[dit()];
    FArrayBox& muUFab = muU[dit()];
    muUFab.mult(muFab, 0, 0);
    dtM = max(muUFab.max(curBox, 0), dtM);
  }
  
  // Gather maximum wave speeds and broadcast the maximum over these
  Vector<Real> allDtMs;
  gather(allDtMs,dtM,uniqueProc(SerialTask::compute));

  if (procID() == uniqueProc(SerialTask::compute)) {
    Real eps = 1e-20;
    dtM = allDtMs[0];
    for (int i = 1; i < allDtMs.size (); ++i)
        dtM = Max(dtM,allDtMs[i]);
    dtM = numerical::Ad * 1/max(dtM, eps);
  }
  
  broadcast(dtM,uniqueProc(SerialTask::compute));
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::computeDtM " << m_level << " dtM " << dtM*tBar << endl;
  return dtM;
}

/*******/
Real
AMRLevelAdvectDiffuse::
computeInitialDt()
{
  Real newDT, eps = 1e-20;
  double speedMax;
  if((speedMax = m_levelGodunov.getMaxWaveSpeed(m_UNew, m_advVel)) > eps)
    newDT = m_initial_dt_multiplier * m_dx / speedMax;
  else if (m_nu > eps)
    newDT = m_cfl*m_dx*m_dx/m_nu;
  else
    MayDay::Error("Insufficient condition to determine time step");
  
  return newDT;
}

/*******/
void
AMRLevelAdvectDiffuse::
levelSetup()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectDiffuse::levelSetup " << m_level << endl;
  }
  
  AMRLevelAdvectDiffuse* amrADCoarserPtr = getCoarserLevel();
  AMRLevelAdvectDiffuse* amrADFinerPtr   = getFinerLevel();
  
  m_hasCoarser = (amrADCoarserPtr != NULL);
  m_hasFiner   = (amrADFinerPtr   != NULL);
  
  int nRefCrse = -1;
  DisjointBoxLayout* crseGridsPtr = NULL;
  
  if (m_hasCoarser)
  {
    nRefCrse = m_coarser_level_ptr->refRatio();
    crseGridsPtr = &amrADCoarserPtr->m_grids;
    
    const DisjointBoxLayout& coarserLevelDomain = amrADCoarserPtr->m_grids;
    
    
    
    m_coarseAverage.define(m_grids,
                           1,
                           nRefCrse);
    
    m_coarseAverageM.define(m_grids,
                           m_gas.m_numOfIonSpe,
                           nRefCrse);
    
    m_coarseAverageFace.define(m_grids,
                           1,
                           nRefCrse);
    
    m_fineInterp.define(m_grids,
                        1,
                        nRefCrse,
                        m_problem_domain);
    
    m_fineInterpM.define(m_grids,
                        m_gas.m_numOfIonSpe,
                        nRefCrse,
                        m_problem_domain);
    
    m_fineInterpVec.define(m_grids,
                           SpaceDim,
                           nRefCrse,
                           m_problem_domain);
    
    m_fineInterpFace.define(m_grids,
                            1,
                            nRefCrse,
                            m_problem_domain);
    
    m_fineInterpPhtzn.define(m_grids,
                           m_phtzn.Psi.nComp(),
                           nRefCrse,
                           m_problem_domain);
    
    m_pwl.define(m_grids, amrADCoarserPtr->m_grids, 1, amrADCoarserPtr->m_problem_domain, amrADCoarserPtr->m_ref_ratio, m_numGhost);
    
    m_pwlM.define(m_grids, amrADCoarserPtr->m_grids, m_gas.m_numOfIonSpe, amrADCoarserPtr->m_problem_domain, amrADCoarserPtr->m_ref_ratio, m_numGhost);
    
    m_pwlVec.define(m_grids, amrADCoarserPtr->m_grids, SpaceDim, amrADCoarserPtr->m_problem_domain, amrADCoarserPtr->m_ref_ratio, m_numGhost);
    
    m_pwlf.define(m_grids, amrADCoarserPtr->m_grids, 1, amrADCoarserPtr->m_problem_domain, amrADCoarserPtr->m_ref_ratio, m_numGhost);
    
    m_qcfi.define(m_grids, &(amrADCoarserPtr->m_grids), m_dx, amrADCoarserPtr->m_ref_ratio, 1, m_problem_domain);
    
    // Maintain levelGodunov
    m_levelGodunov.define(*m_advPhys,
                          m_grids,
                          coarserLevelDomain,
                          m_problem_domain,
                          nRefCrse,
                          m_useLimiting,
                          m_dx,
                          m_hasCoarser,
                          m_hasFiner);
    
    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.
    
    // Maintain flux registers
    amrADCoarserPtr->m_fluxRegister.define(m_grids,
                                           amrADCoarserPtr->m_grids,
                                           m_problem_domain,
                                           amrADCoarserPtr->m_ref_ratio,
                                           1);
    amrADCoarserPtr->m_fluxRegister.setToZero();

  }
  else
  {
    m_levelGodunov.define(*m_advPhys,
                          m_grids,
                          DisjointBoxLayout(),
                          m_problem_domain,
                          m_ref_ratio,
                          m_useLimiting,
                          m_dx,
                          m_hasCoarser,
                          m_hasFiner);
  }
  
  defineSolvers();
}

/*******/
AMRLevelAdvectDiffuse*
AMRLevelAdvectDiffuse::
getCoarserLevel() const
{
  AMRLevelAdvectDiffuse* amrADCoarserPtr = NULL;
  
  if (m_coarser_level_ptr != NULL)
  {
    amrADCoarserPtr = dynamic_cast<AMRLevelAdvectDiffuse*>(m_coarser_level_ptr);
    
    if (amrADCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelAdvectDiffuse::getCoarserLevel: dynamic cast failed");
    }
  }
  
  return amrADCoarserPtr;
}

/*******/
AMRLevelAdvectDiffuse*
AMRLevelAdvectDiffuse::
getFinerLevel() const
{
  AMRLevelAdvectDiffuse* amrADFinerPtr = NULL;
  
  if (m_finer_level_ptr != NULL)
  {
    amrADFinerPtr = dynamic_cast<AMRLevelAdvectDiffuse*>(m_finer_level_ptr);
    
    if (amrADFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelAdvectDiffuse::getFinerLevel: dynamic cast failed");
    }
  }
  
  return amrADFinerPtr;
}

/*******/
void
AMRLevelAdvectDiffuse::
fillMobility(bool timeInterpForGhost) {
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::fillMobility " << m_level << endl;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
    FArrayBox& uFab = m_mu[dit()];
    FArrayBox& EmagFab = m_field.m_Emag[dit()];
    // don't deal with ghost cells
    const Box& b = m_grids[dit()];
//    const Box& fabBox = uFab.box();
    for (BoxIterator bit(b); bit.ok(); ++bit) {
      const IntVect& iv = bit();
      Real n;
      if (m_gas.m_uniformity)
        n = m_gas.m_N;
      else
        n = m_neut[dit()](iv,0);
      
      uFab(iv, 0) = m_gas.EDpdentProcs["mobility"].value(EmagFab(iv,0)/n)/n;
    }
  }
  
  Real eps = 1.0e-10;
  AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
  if (m_hasCoarser) {
    if (!timeInterpForGhost)
      fillGhost(m_pwl, m_mu, amrGodCoarserPtr->m_mu, 0, amrGodCoarserPtr->m_mu, 0, eps, 0);
    else
      fillGhost(m_pwl, m_mu, amrGodCoarserPtr->m_muOld, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_mu, amrGodCoarserPtr->m_time, m_dt, m_time+m_dt);
  } else
    m_mu.exchange();
  
  averageCellToEdge(m_grids, m_mu, m_muEdge);
    
  if (s_verbosity >= 3) {
    printDiagnosticInfo (m_level, m_dx, m_grids, m_mu, "mob", "AMRLevelAdvectDiffuse::fillMobility ");
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
fillAdvectionVelocity(bool timeInterpForGhost) {
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::fillAdvectionVelocity " << m_level << endl;
  
//  initialization
//  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
//    m_advVel[dit()].setVal(0.0, m_advVel[dit()].box());
  
  m_field.m_EEdge.copyTo(m_advVel);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
//    Box b = m_grids[dit()];
    Box b = m_advVel[dit()].box();
    m_advVel[dit()].mult(m_muEdge[dit()], b, 0, 0);
    m_advVel[dit()].negate();
    pout() << "advVelBox: " << b << " muEdgBox: " << m_muEdge[dit()].box() << endl;
  }
  
  Real eps = 1.0e-10;
  AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
  if (m_hasCoarser) {
    if (!timeInterpForGhost)
      fillGhostFace(m_advVel, amrGodCoarserPtr->m_advVel, 0, amrGodCoarserPtr->m_advVel, 0, eps, 0);
    else
      fillGhostFace(m_advVel, amrGodCoarserPtr->m_advVelOld, amrGodCoarserPtr->m_time-amrGodCoarserPtr->m_dt, amrGodCoarserPtr->m_advVel, amrGodCoarserPtr->m_time, m_dt, m_time+m_dt);
  } else
    m_advVel.exchange();
  
  if (s_verbosity >= 3) {
    printDiagnosticInfo (m_level, m_dx, m_grids, m_advVel, "advVel", "fillAdvectionVelocity");
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
testing () {
  if (s_verbosity >= 3)
    pout() << "AMRLevelAdvectDiffuse::testing " << m_level << endl;
  
//  // test getDivDeltaFlux
//  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
//    m_testOutput[dit()].setVal(0.0, m_testOutput[dit()].box(), 0);
//    m_field.m_EEdge[dit()].setVal(0.0, m_field.m_EEdge[dit()].box());
//    m_field.m_EEdge[dit()].setVal((m_level+1)*(m_level+1));
//  }
//
//  if (m_hasFiner) {
//    AMRLevelAdvectDiffuse* fl = getFinerLevel();
//    getDivDeltaFlux(m_testOutput, &(m_field.m_EEdge), &(fl->m_field.m_EEdge), NULL, 1.0);
//  }
  
  if (m_gas.m_uniformity) {
    getRate(m_testOutput, m_field.m_Emag, m_mu, "ionization");
  } else {
    getRate(m_testOutput, m_field.m_Emag, m_neut, m_mu, "ionization");
  }
}

/*******/
void
AMRLevelAdvectDiffuse::
outputStepStats(std::ofstream& ofs) {
  
  unsigned long long localNumAdvCells = 0;
  unsigned long long currNumAdvCells = 0;
  Real time_eps = 1.0e-20;
  
  if (m_time > time_eps) {
    localNumAdvCells = m_grids.numPointsThisProc();
    // cout << "cpu = " << procID() << " level = " << m_level << " localNumAdvCells = " << m_grids.numPointsThisProc() << endl;
    
    // Gather and broadcast
    Vector<unsigned long long> allLocalNumAdvCells;
    gather(allLocalNumAdvCells,localNumAdvCells,uniqueProc(SerialTask::compute));

    if (procID() == uniqueProc(SerialTask::compute)) {
      for(int i = 0; i < allLocalNumAdvCells.size(); ++i)
        currNumAdvCells += allLocalNumAdvCells[i];
      totNumAdvCells += currNumAdvCells;
    }
    
    broadcast(totNumAdvCells,uniqueProc(SerialTask::compute));
    broadcast(currNumAdvCells,uniqueProc(SerialTask::compute));
    // cout << "AMRLevelAdvectDiffuse::outputStepStats " << m_level << " totNumAdvCells = " << totNumAdvCells << endl;
  }
  
  Vector<AMRLevelAdvectDiffuse*>         hierarchy;
  Vector<int>                            refRat;
  Vector<DisjointBoxLayout>              grids;
  Real                                   lev0Dx;
  ProblemDomain                          lev0Domain;
  getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);
  int finest_level = hierarchy.size()-1;
  
  unsigned long long totNumCells = 0;
  unsigned long long localNumPts = 0;
  for (int lev = 0; lev <= finest_level; lev++) {
    localNumPts += grids[lev].numPointsThisProc();
  //  cout << "cpu = " << procID() << " level = " << lev << " time = " << hierarchy[lev]->time() << " localLevelNumPts = " << grids[lev].numPointsThisProc() << endl;
  }
  
  // Gather and broadcast
  Vector<unsigned long long> allLocalNumPts;
  gather(allLocalNumPts,localNumPts,uniqueProc(SerialTask::compute));

  if (procID() == uniqueProc(SerialTask::compute)) {
    for(int i = 0; i < allLocalNumPts.size(); ++i)
         totNumCells += allLocalNumPts[i];
  }
  
  broadcast(totNumCells,uniqueProc(SerialTask::compute));
  
  if (procID() == uniqueProc(SerialTask::compute)) {
    
    if (m_time < time_eps && m_level == finest_level) {
      ofs << setiosflags(ios::left) << setw(12) << "cpuTime";
      ofs << setiosflags(ios::left) << setw(6) << "lev";
      ofs << setiosflags(ios::left) << setw(12) << "cells";
      ofs << setiosflags(ios::left) << setw(12) << "advCells";
      ofs << setiosflags(ios::left) << setw(12) << "chgAdvCells";
      ofs << setiosflags(ios::left) << setw(12) << "levTimes";
      ofs << endl;
    }
#ifdef CH_MPI
    ofs << setiosflags(ios::left) << setw(12) << MPI_Wtime() - startWTime;
#endif
//    ofs << setiosflags(ios::left) << setw(12) << timer.getTimeStampWC();
    ofs << setiosflags(ios::left) << setw(6) << m_level;
    ofs << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << totNumCells;
    ofs << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << totNumAdvCells;
    ofs << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << currNumAdvCells;
    for (int lev = 0; lev <= finest_level; lev++) {
      ofs << setiosflags(ios::left) << setw(12) << setprecision(4) << hierarchy[lev]->time();
    }
    ofs << endl;
    
  }
}


#include "NamespaceFooter.H"
