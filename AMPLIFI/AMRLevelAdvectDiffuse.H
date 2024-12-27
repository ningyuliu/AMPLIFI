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

#ifndef _AMRLEVELADVECTDIFFUSE_H_
#define _AMRLEVELADVECTDIFFUSE_H_
#include <fstream>

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "FineInterp.H"
#include "FineInterpFace.H"
#include "LevelFluxRegister.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "QuadCFInterp.H"

#include "LevelAdvect.H"
#include "PhysIBC.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "RealVect.H"
#include "LevelFluxRegister.H"
#include "DisjointBoxLayout.H"
#include "AdvectionFunctions.H"
#include "BCFunc.H"
#include "LevelTGA.H"
#include "AMRPoissonOp.H"
#include "BiCGStabSolver.H"
#include "RelaxSolver.H"

#include "physicalConstants.h"
#include "gas.hpp"
#include "electricField.hpp"
#include "photoionization.hpp"

#include "NamespaceHeader.H"

/// AMRLevel for advection-diffusion
/**
 */
class AMRLevelAdvectDiffuse : public AMRLevel
{
public:
  
  /// Default constructor
  AMRLevelAdvectDiffuse()
  {
    m_isDefined = false;
  }
  
  /// Full constructor. Arguments are same as in define()
  AMRLevelAdvectDiffuse(const AdvectPhysics&        a_gphys,
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
    define(a_gphys, a_gas, a_bcFunc, a_EPotbcFunc, a_PIbcFunc, a_cfl, a_domainLength, a_refineThresh,
           a_tagBufferSize, a_initialDtMultiplier, a_useLimiting, a_nu);
  }
    
  void getRate(LevelData<FArrayBox>& a_rate, const LevelData<FArrayBox>& a_Emag, const LevelData<FArrayBox>& a_mu, string procName);
  void getRate(FArrayBox& a_rate, const FArrayBox& a_Emag, const FArrayBox& a_mu, string procName);
  
  void getRate(LevelData<FArrayBox>& a_rate, const LevelData<FArrayBox>& a_Emag, const LevelData<FArrayBox>& a_n, const LevelData<FArrayBox>& a_mu, string procName);
  void getRate(FArrayBox& a_rate, const FArrayBox& a_Emag, const FArrayBox& a_n, const FArrayBox& a_mu, string procName);
  void updateWithReactionContribution (LevelData<FArrayBox>& U, const LevelData<FArrayBox>& Emag, const LevelData<FArrayBox>& mu, const double dt);
  void updateWithReactionContribution (LevelData<FArrayBox>& U, const LevelData<FArrayBox>& Emag, const LevelData<FArrayBox>& mu, const LevelData<FArrayBox>& phtznRate, const double dt);
  void makeDiffusiveSource(LevelData<FArrayBox>& a_diffusiveSource, const LevelData<FArrayBox>& U);
  
  /// Defines this AMRLevelAdvectDiffuse
  void define(/// advection physics class
              const AdvectPhysics&        a_gphys,
              gas                         a_gas,
              /// boundary condition class for diffusion solve
              BCHolder                    a_bcFunc,
              /// boundary condition class for electric field solve
              BCHolder                    a_EPotbcFunc,
              /// boundary condition class for electric field solve
              BCHolder                    a_PIbcFunc,
              /// CFL number
              const Real&                 a_cfl,
              /// physical length of domain
              const Real&                 a_domainLength,
              /// undivided gradient size over which a cell will be tagged for refinement
              const Real&                 a_refineThresh,
              /// number of buffer cells around each tagged cell that will also be tagged
              const int&                  a_tagBufferSize,
              /// CFL number at beginning of calculation
              const Real&                 a_initialDtMultiplier,
              /// whether to use van Leer limiting
              const bool&                 a_useLimiting,
              /// diffusion coefficient
              const Real&                 a_nu);
  
  ///
  virtual ~AMRLevelAdvectDiffuse();
  
  /// Never called: historical
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio)
  {
    MayDay::Error("never called--historical");
  }
  
  /// Define new AMRLevelAdvectDiffuse from coarser
  /**
   */
  virtual void define(
                      AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);
  
  /// Advance by one timestep
  virtual Real advance();
  
  /// Things to do after a timestep
  virtual void postTimeStep();
  
  /// Create tags for regridding
  virtual void tagCells(IntVectSet& a_tags) ;
  
  /// Create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;
  
  /// Set up data on this level after regridding
  virtual void regrid(const Vector<Box>& a_newGrids);
  
  /// Initialize grids
  virtual void initialGrid(const Vector<Box>& a_newGrids);
  
  /// Initialize data
  virtual void initialData();
  
  /// Things to do after initialization
  virtual void postInitialize();
  
#ifdef CH_USE_HDF5
  /// Write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;
  
  /// Write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;
  
  /// Read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);
  
  /// Read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);
  
  /// Write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;
  
  /// Write plotfile data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif
  
  /// Returns the dt computed earlier for this level
  virtual Real computeDt();
  
  virtual Real computeDtI();
  
  virtual Real computeDtM();
  
  /// Compute dt using initial data
  virtual Real computeInitialDt();
  
  //for convergence tests
  LevelData<FArrayBox>& getStateNew()
  {
    return m_UNew;
  }
  
  //for convergence tests
  LevelData<FArrayBox>& getStateOld()
  {
    return m_UOld;
  }
protected:
  void setSolverCoef(Real a_alpha, Real a_beta);
  void interpolateInTime(LevelData<FArrayBox>&          a_interp,
                         const LevelData<FArrayBox>&    a_old,
                         const LevelData<FArrayBox>&    a_new,
                         Real a_time, Real a_tOld, Real a_tNew);
  
  void getHierarchyAndGrids(Vector<AMRLevelAdvectDiffuse*>&        a_hierarchy,
                            Vector<DisjointBoxLayout>&             a_grids,
                            Vector<int>&                           a_refRat,
                            ProblemDomain&                         a_lev0Dom,
                            Real&                                  a_lev0Dx);
  
  void doImplicitReflux();
  
  void printRefluxRHSMax(const std::string& a_preflix) ;
  
  Real diffusiveAdvance(LevelData<FArrayBox>& a_diffusiveSource);
  Real strangAdvance(LevelData<FArrayBox>& a_diffusiveSource);
  
  void defineSolvers();
  BCHolder                    m_bcFunc;
  // solver to do level diffusion and flux register interaction
  static RefCountedPtr<LevelTGA>                                  s_diffuseLevTGA;
  // one diffusionSolver to rule them all and in darkness bind them
  static RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      s_diffuseAMRMG;
  // operator factory
  static RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > s_diffuseOpFact;
  // bottom solver
  static BiCGStabSolver<LevelData<FArrayBox> >                    s_botSolver;
  // static RelaxSolver<LevelData<FArrayBox> >                       s_botSolver;
  
  BCHolder                    m_EPotbcFunc;
  // one AMRMG solver to find electric potential for the explicit scheme
  static RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      s_EPotAMRMG;
  // operator factory
  static RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > s_EPotOpFact;
  // bottom solver
  static BiCGStabSolver<LevelData<FArrayBox> >                    s_EPotBotSolver;
  
  // there was an issue in using one AMRMG solver to find electric potential
  // for the implicit scheme. So the solver is generated on the fly at the solve
  // step.
//  static RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      s_EPotImpAMRMG;
//  // operator factory
//  static RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > s_EPotImpOpFact;
  // bottom solver
  static BiCGStabSolver<LevelData<FArrayBox> >                    s_EPotImpBotSolver;
  static BiCGStabSolver<LevelData<FArrayBox> >                    s_EPotImpCompBotSolver;
  
  BCHolder                    m_phtznbcFunc;
  // one AMRMG solver for photoionization
  static RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      s_PIAMRMG;
  // operator factory
  static RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > s_PIOpFact;
  // bottom solver
  static BiCGStabSolver<LevelData<FArrayBox> >                    s_PIBotSolver;
  
  //embarrassing but needed function
  void getCoarseDataPointers(LevelData<FArrayBox>** a_coarserDataOldPtr,
                             LevelData<FArrayBox>** a_coarserDataNewPtr,
                             LevelFluxRegister**    a_coarserFRPtr,
                             LevelFluxRegister**    a_finerFRPtr,
                             Real& a_tCoarserOld,
                             Real& a_tCoarserNew);
  
  void fillGhost(PiecewiseLinearFillPatch&    a_pwl,
                 LevelData<FArrayBox>&       a_U,
                 const LevelData<FArrayBox>& a_UCoarseOld,
                 const Real&                 a_TCoarseOld,
                 const LevelData<FArrayBox>& a_UCoarseNew,
                 const Real&                 a_TCoarseNew,
                 const Real&                 a_dt,
                 const Real&                 a_time);
  
  void fillGhostFace(LevelData<FluxBox>&       a_U,
                     const LevelData<FluxBox>& a_UCoarseOld,
                     const Real&                 a_TCoarseOld,
                     const LevelData<FluxBox>& a_UCoarseNew,
                     const Real&                 a_TCoarseNew,
                     const Real&                 a_dt,
                     const Real&                 a_time);
  
  void fillGhostQuad(QuadCFInterp&               a_qcfi,
                     LevelData<FArrayBox>&       a_U,
                     const LevelData<FArrayBox>& a_UCoarseOld,
                     const Real&                 a_TCoarseOld,
                     const LevelData<FArrayBox>& a_UCoarseNew,
                     const Real&                 a_TCoarseNew,
                     const Real&                 a_time);
  
  // solve Poisson's equation on signle level
  void poissonSolve();
  // solve Poisson's equation across all levels
  void poissonSolveComposite();
  // solve Poisson's equation with variable coefficients on signle level
  void poissonSolveImplicit();
  // for variable coeff Poisson:
  // (alpha * aCoef(x) * I - beta * Div(bCoef(x) . Grad)) phi = rho
  // this sets aCoef and bCoef
  void setPoissonCoeffAB();
  // solve Poisson's equation with variable coefficients across all levels
  void poissonSolveImplicitComposite();
  void setPoissonCoeffABComposite(Real commonDt);
  // compute div(delta F) see Chombo Design document or Colella et al. 1999a,b
  void getDivDeltaFlux(LevelData<FArrayBox>& rhs, LevelData<FluxBox>* edgeField, LevelData<FluxBox>* edgeFieldFine, LevelData<FluxBox>* coeff, Real scale = 1.0);
  void computeEField(bool timeInterpForGhost);
  void photoionizationSolve();
  
  void fillMobility(bool timeInterpForGhost);
  void fillAdvectionVelocity(bool timeInterpForGhost);
  // Setup menagerie of data structures
  void levelSetup();
  
  // Get the next coarser level
  AMRLevelAdvectDiffuse* getCoarserLevel() const;
  
  // Get the next finer level
  AMRLevelAdvectDiffuse* getFinerLevel() const;
  
  void testing();
  void outputStepStats (std::ofstream&);
  
  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld, m_UNew, m_dU, m_ionNew, m_ionOld, m_neut;
  LevelData<FluxBox>   m_advVel, m_advVelOld;
  
  Vector<string> m_stateNames;
  bool m_isDefined;
  Real                 m_cfl;
  Real                 m_domainLength;
  Real                 m_refineThresh;
  int                  m_tagBufferSize;
  Real                 m_initialDtMultiplier;
  bool                 m_useLimiting;
  Real                 m_nu;
  Real                 m_dx;
  RefCountedPtr<AdvectPhysics> m_advPhys;
  AdvectionVelocityFunction    m_advFunc;
  
  FineInterp m_fineInterp;
  FineInterp m_fineInterpM;
  FineInterp m_fineInterpVec;
  FineInterpFace m_fineInterpFace;
  CoarseAverage m_coarseAverage;
  CoarseAverage m_coarseAverageM;
  CoarseAverageFace m_coarseAverageFace;
  PiecewiseLinearFillPatch m_pwl;
  PiecewiseLinearFillPatch m_pwlM;
  PiecewiseLinearFillPatch m_pwlVec;
  PiecewiseLinearFillPatchFace m_pwlf;
  QuadCFInterp m_qcfi;
  
  Real m_dtNew;
  int m_numGhost;
  LevelAdvect m_levelGodunov;
  LevelFluxRegister m_fluxRegister;
  
  bool m_hasCoarser;
  bool m_hasFiner;
  bool m_doImplicitReflux;
  bool m_hasDiffusion;
  
  DisjointBoxLayout m_grids;
  
  gas m_gas;
  // m_fieldOld and m_phiOld may not be needed;
  field m_field, m_fieldOld;
  LevelData<FArrayBox> m_phiOld, m_phi;
  LevelData<FArrayBox> m_mu;
  LevelData<FArrayBox> m_muOld;
  photoionization m_phtzn, m_phtznOld;
  FineInterp m_fineInterpPhtzn;

  LevelData<FArrayBox> m_dUDiff;
  LevelData<FluxBox>   m_flux, m_muEdge;
  
  bool m_varyingField, m_doImplicitPoisson;
  static Vector<RefCountedPtr<LevelData<FArrayBox> > > s_aCoef;
  static Vector<RefCountedPtr<LevelData<FluxBox> > > s_bCoef;
  static Vector<RefCountedPtr<LevelData<FArrayBox> > > s_aCoefComp;
  static Vector<RefCountedPtr<LevelData<FluxBox> > > s_bCoefComp;
  LevelData<FluxBox>   m_EEdgeOld, m_EEdgeCoarseOld;
  
  int m_sourceNumericalScheme;
  
  static bool s_testing;
  LevelData<FArrayBox> m_testOutput;
  static unsigned long long totNumAdvCells;
private:
  
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelAdvectDiffuse&);
  AMRLevelAdvectDiffuse(const AMRLevelAdvectDiffuse&);
};

#include "NamespaceFooter.H"

#endif
