//
//  photoionization.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef photoionization_hpp
#define photoionization_hpp

#include "RealVect.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"

#include "NamespaceHeader.H"

#include "physicalConstants.h"
#include "globalVariables.h"

using namespace normalization;

// photoionization (PI) SP3 model rates and coefficients
//const double pO2Torr = 150;
//const double PIA1 = 6.7e-3*1e2*pO2Torr*lBar; // 1/(m Torr);
//const double PIA2 = 0.0346*1e2*pO2Torr*lBar; // 1/(m Torr);
//const double PIA3 = 0.3059*1e2*pO2Torr*lBar; // 1/(m Torr);
//const double PIlambda1 = 0.0447*1e2*pO2Torr*lBar; // 1/m/Torr
//const double PIlambda2 = 0.1121*1e2*pO2Torr*lBar; // 1/m/Torr
//const double PIlambda3 = 0.5994*1e2*pO2Torr*lBar; // 1/m/Torr
const double pO2Torr = 150;
const double quenchingFactor = 30.0/(760.0+30.0);
const double PIXi = 0.06;
extern vector<double> SP3A;
extern vector<double> SP3Lambda;

//  -----------------------------------------
// boundary condition stuff
//  -----------------------------------------
// use the same integer array to store the types of the BCs for all 3 SP3 components
extern vector<int>    PIBCLo;
extern vector<int>    PIBCHi;
extern double         PIDiriNeumBCValues;
// Robin BC: (alpha + d/dx) phi = gamma
// store alphaLo, alphaHi, and gamma in order for a SP3 component;
// also they are assumed to be the same for all directions
extern vector<double> PIRobinBCValues;

// get and set the parameters of the multigrid solver for PI
void PISetupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver);

// functions for testing the PI model (see Bourdon et al. (2007)).
enum class PITestProbType {zeroRHS = 0, unityRHS, gaussian};
void PISetTestRHS(Vector<LevelData<FArrayBox>* > a_rhs, Vector<int>& a_refRat, ProblemDomain a_lev0Dom, Real a_lev0Dx, int startLev, int a_finestLevel, PITestProbType s_probtype);

void PIParseValue(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);
void PIParseValueRobin(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);
void PIParseBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous);

class photoionization {
  
public:
  photoionization() {A.reserve(10); lambda.reserve(10);}
  void define(DisjointBoxLayout a_grids, int ncomps, IntVect ivGhost);
  void define(const photoionization& a_pi);
  void copyTo(photoionization& des);
  void setCoefficients(vector<double> a_A, vector<double> a_lambda);
  void calcRate();
  virtual ~photoionization() {};
  
  bool runSolve;
  vector<double> A;
  vector<double> lambda;
  int ncomps;
  LevelData<FArrayBox> Psi, rate;
 
protected:
  // verbosity level
  static int s_verbosity;
  
};

#include "NamespaceFooter.H"

#endif /* photoionization_hpp */
