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

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <functional>
#include <cfloat>

#include "AdvectDiffuseUtils.H"
#include "ParmParse.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "CellToEdge.H"
#include "electricField.hpp"
#include "globalVariables.h"
#include "parameterizedFunction.hpp"
#include "electricField.hpp"
#include "blob.hpp"

#include "NamespaceHeader.H"

// give default values
namespace normalization {

double scalingFactor = 1.0;

double EBar = 1e5;
double muBar = 0.0382;
double lBar = sqrt(constants::e/(constants::eps0*EBar));
double tBar = lBar/(muBar*EBar);
double nBar = pow(lBar, -3);
double phiBar = EBar*lBar;

}

namespace numerical {

double Ac;
double Ai;
double Ad;

}

std::string       AMPLIFIOutFilename("AMPLIFIOutput.txt");
std::ofstream     AMPLIFIOut;
double            startWTime;
OldTimer          timer;

std::string       timeSeriesFilename("timeSeries.txt");
std::ofstream     timeSeriesOut;

void ADParseValue(Real* pos,
                  int* dir,
                  Side::LoHiSide* side,
                  Real* a_values)
{
  ParmParse pp("diffusion");
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}

void ADParseBC(FArrayBox& a_state,
               const Box& a_valid,
               const ProblemDomain& a_domain,
               Real a_dx,
               bool a_homogeneous) {
  if (!a_domain.domainBox().contains(a_state.box())) {
    
    Box valid = a_valid;
    int a_order = 1;
    
    for (int i=0; i<CH_SPACEDIM; ++i) {
      // don't do anything if periodic
      if (!a_domain.isPeriodic(i)) {
        ParmParse pp("diffusion");
        std::vector<int>  bcLo = std::vector<int>();
        std::vector<int>  bcHi = std::vector<int>();
        pp.getarr("bc_lo", bcLo, 0, SpaceDim);
        pp.getarr("bc_hi", bcHi, 0, SpaceDim);
        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        
        if (!a_domain.domainBox().contains(ghostBoxLo))
          switch (bcLo[i]) {
            case 0:
              DiriBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     ADParseValue,
                     i,
                     Side::Lo);
              break;
              
            case 1:
              NeumBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     ADParseValue,
                     i,
                     Side::Lo);
              break;
              
            case 2:
              ExtrapolateBC(a_state,
                            valid,
                            a_dx,
                            i,
                            Side::Lo,
                            a_order);
              break;
              
            default:
              MayDay::Error("bogus bc flag lo");
              break;
          }
        
        if (!a_domain.domainBox().contains(ghostBoxHi))
          switch (bcHi[i]) {
            case 0:
              DiriBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     ADParseValue,
                     i,
                     Side::Hi);
              break;
              
            case 1:
              NeumBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     ADParseValue,
                     i,
                     Side::Hi);
              break;
              
            case 2:
              ExtrapolateBC(a_state,
                            valid,
                            a_dx,
                            i,
                            Side::Hi,
                            a_order);
              break;
              
            default:
              MayDay::Error("bogus bc flag ho");
              break;
          }
      } // end if is not periodic in ith direction
    }
  }
}

int
orderScript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}

void
getErrorFromCoarseAndFine(Vector< LevelData<FArrayBox>* >&           a_errorCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnCoar,
                          const Vector< DisjointBoxLayout >&         a_gridsCoar,
                          const ProblemDomain&                       a_level0DomainCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnFine,
                          const Vector< DisjointBoxLayout >&         a_gridsFine,
                          const ProblemDomain&                       a_level0DomainFine,
                          const Vector<int>&                         a_refRat)
{
  
  int nlevels = a_solnFine.size();
  a_errorCoar.resize(nlevels);
  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two solutions
  int nvar = a_solnFine[0]->nComp();
  Interval interv(0, nvar-1);
  
  ProblemDomain domLevCoar = a_level0DomainCoar;
  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    int nvar = a_solnCoar[ilev]->nComp();
    
    a_errorCoar[ilev] = new LevelData<FArrayBox>(a_gridsCoar[ilev], nvar,  IntVect::Zero);
    
    CoarseAverage averageOp(a_gridsFine[ilev], nvar, nref);
    //here make error = Ave(fine)
    averageOp.averageToCoarse(*a_errorCoar[ilev], *a_solnFine[ilev]);
    //now subtract off coarse so error= Ave(Fine) - coar
    for (DataIterator dit = a_gridsCoar[ilev].dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       errorFAB = (*a_errorCoar[ilev])[dit()];
      const FArrayBox& solnFAB  = (*a_solnCoar[ilev])[dit()];
      
      errorFAB -= solnFAB;
    }
    domLevCoar.refine(a_refRat[ilev]);
  }
}
Real
dtgNorm(const Vector< LevelData<FArrayBox>* >& a_src,
        const Vector< DisjointBoxLayout >&     a_grids,
        const Vector<int>&                     a_refRatio,
        const ProblemDomain&                   a_domain,
        const int& a_comp,
        const int& a_pval)
{
  Real dx = 1;
  Interval interv (a_comp, a_comp);
  int lbase = 0;
  Real volWeightedSum = computeNorm(a_src, a_refRatio, dx, interv, a_pval, lbase);
  
  //now unweight it from the volume.
  Real norm = volWeightedSum;
  if (a_pval != 0)
  {
    Real volume = 1;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      volume *= a_domain.size(idir);
    }
    norm /= volume;
  }
  
  return norm;
}
/******/
void
compareError(Vector<Real>&                            a_orders,
             const Vector< LevelData<FArrayBox>* >&   a_errorFine,
             const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
             const Vector< DisjointBoxLayout >&       a_gridsFine,
             const Vector< DisjointBoxLayout >&       a_gridsCoar,
             const Vector<int>&                       a_refRat,
             const ProblemDomain &                    a_coarseDom,
             int a_testverbosity)
{
  ProblemDomain fineDom = a_coarseDom;
  fineDom.refine(2);
  const Vector<int>& refRat = a_refRat;
  const int ncomp = a_errorFine[0]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  a_orders.resize(ncomp*nnorm, 0.0);
  Real* orders    = &(a_orders[0]);
  for (int icomp = 0; icomp < ncomp; icomp++)
  {
    orders[icomp] = 0.0;
    for (int inorm = 0; inorm < nnorm; inorm++)
    {
      normsCoar[orderScript(icomp, inorm, ncomp)] = 0;
      normsFine[orderScript(icomp, inorm, ncomp)] = 0;
    }
  }
  int testverbosity = a_testverbosity;
  for (int comp = 0; comp < ncomp; comp++)
  {
    for (int inorm = 0; inorm <= 2; inorm++)
    {
      
      Real coarnorm = dtgNorm(a_errorCoar,
                              a_gridsCoar,
                              refRat, a_coarseDom,
                              comp, inorm);
      Real finenorm = dtgNorm(a_errorFine,
                              a_gridsFine,
                              refRat, fineDom,
                              comp, inorm);
      normsCoar[orderScript(comp,inorm,ncomp)] = coarnorm;
      normsFine[orderScript(comp,inorm,ncomp)] = finenorm;
      Real order = 0;
      if ((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
      {
        order = log(Abs(coarnorm/finenorm))/log(2.0);
      }
      orders[orderScript(comp,inorm,ncomp)] = order;
    }
  }
  
  //output in latex format
  int nfine = a_coarseDom.size(0);
  nfine *= 2;
  if (testverbosity > 0)
  {
    pout() << setw(12)
    << setprecision(6)
    << setiosflags(ios::showpoint)
    << setiosflags(ios::scientific) ;
    for (int inorm = 0; inorm <= 2; inorm++)
    {
      pout() << "\\begin{table}[p]" << endl;
      pout() << "\\begin{center}" << endl;
      pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
      pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
      pout() << "\\hline \\hline " << endl;
      for (int icomp = 0; icomp < ncomp; icomp++)
      {
        int iindex = orderScript(icomp,inorm,ncomp);
        pout() << "var" << icomp << " &    \t "
        << setw(12)
        << setprecision(6)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << normsCoar[iindex]  << " & "
        << setw(12)
        << setprecision(6)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << normsFine[iindex] << " & "
        << setw(12)
        << setprecision(6)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << orders[iindex];
        pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
      }
      pout() << "\\end{tabular}" << endl;
      pout() << "\\end{center}" << endl;
      pout() << "\\caption{Solution error convergence rates using L-" << inorm << " norm. " << endl;
      pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
      pout() << "\\end{table}" << endl;
      pout() << endl << endl;
    }
  }
  //latex output
  
  delete[] normsFine;
  delete[] normsCoar;
}
void
makeFinestDomain(ProblemDomain& a_domain,
                 Real&          a_dx)
{
  ParmParse pp;
  Real domainLength;
  pp.get("domain_length",domainLength);
  domainLength = domainLength/normalization::scalingFactor/normalization::lBar;
  getProblemDomain(a_domain);
  int ncells = a_domain.size(0);
  a_dx = domainLength/ncells;
}
///
/**
 **/
void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar)
{
  a_boxesCoar.resize(a_boxesFine.size());
  for (int ilev = 0; ilev < a_boxesFine.size(); ilev++)
  {
    a_boxesCoar[ilev].resize(a_boxesFine[ilev].size());
    for (int ibox = 0; ibox < a_boxesFine[ilev].size(); ibox++)
    {
      a_boxesCoar[ilev][ibox] = coarsen(a_boxesFine[ilev][ibox], a_refToCoar);
    }
  }
}
/*****/
void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain)
{
  ParmParse pp;
  int maxLevel;
  pp.get("max_level", maxLevel);
  if (maxLevel == 1)
  {
    int amrRef = 2;
    a_refRat.resize(2, amrRef);
    a_boxes.resize(2);
    Box fineBox = refine(a_domain, amrRef);
    int ishrink = fineBox.size(0);
    //this leaves 1/4 refined.
    ishrink *= 3;
    ishrink /= 8;
    fineBox.grow(-ishrink);
    a_boxes[0] = Vector<Box>(1, a_domain);
    a_boxes[1] = Vector<Box>(1, fineBox);
  }
  else if (maxLevel == 0)
  {
    a_refRat.resize(1, 2);
    a_boxes.resize(1);
    a_boxes[0] = Vector<Box>(1, a_domain);
  }
  else
  {
    MayDay::Error("can only deal with two levels now");
  }
}

void
getProblemDomain(ProblemDomain& a_domain)
{
  ParmParse pp;
  Vector<int> ncell(SpaceDim);
  pp.getarr("n_cell", ncell, 0, SpaceDim);
  IntVect hiEnd(D_DECL(ncell[0]-1,ncell[1]-1,ncell[2]-1));
  Box level0Domain(IntVect::Zero, hiEnd);
  
  Vector<int> v_is_periodic(SpaceDim);
  pp.getarr("periodic_bc", v_is_periodic, 0, SpaceDim);
  
  
  bool is_periodic[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++) is_periodic[idir] = (v_is_periodic[idir]==1);
  
  ProblemDomain prob_domain(level0Domain.smallEnd(),
                            level0Domain.bigEnd(),
                            &(is_periodic[0]));
  
  a_domain = prob_domain;
}

void getGlobalVariables() {
  
  ParmParse pp;
  
  Real Ac, Ai, Ad;
  
  pp.get("Ac", Ac);
  pp.get("Ai", Ai);
  pp.get("Ad", Ad);
  
  numerical::Ac = Ac;
  numerical::Ai = Ai;
  numerical::Ad = Ad;
}

void
getAmbientGas(gas& a_gas)
{
  std::string name;
  bool uniformity;
  int numOfIonSpe;
  
  ParmParse pp("gas");
  pp.get("name", name);
  pp.get("numOfIonSpe", numOfIonSpe);
  pp.get("uniformity", uniformity);
  
  Real N, diffCoef;
  pp.get("density", N);
  pp.get("elecDiffCoef", diffCoef);
  
  N *= normalization::scalingFactor;
  N /= normalization::nBar;
  diffCoef /= normalization::scalingFactor;
  diffCoef /= normalization::lBar*normalization::lBar/normalization::tBar;
  
  gas air(name, uniformity, N, numOfIonSpe, diffCoef);
  
  pO2Torr = pO2Torr*normalization::scalingFactor;
  std::transform(SP3A.begin(), SP3A.end(), SP3A.begin(),
                 std::bind(std::multiplies<double>(), std::placeholders::_1, pO2Torr*lBar));
  std::transform(SP3Lambda.begin(), SP3Lambda.end(), SP3Lambda.begin(),
                 std::bind(std::multiplies<double>(), std::placeholders::_1, pO2Torr*lBar));
  
  if (!uniformity) {
    if (pp.contains("densityInputFile")) { // input data from a file
      string inputFileNamePtr;
      pp.get("densityInputFile", inputFileNamePtr);
      Vector<int> numGridPoints(SpaceDim, 0);
      vector<double> spacing(SpaceDim, 0);
      vector<double> basePt(SpaceDim, 0);
      pp.getarr("numGridPoints", numGridPoints, 0, SpaceDim);
      pp.getarr("spacing", spacing, 0, SpaceDim);
      pp.getarr("basePoint", basePt, 0, SpaceDim);
      //non-dimensionalize
      std::transform(spacing.begin(), spacing.end(), spacing.begin(),
                     std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0/lBar));
      //non-dimensionalize
      std::transform(basePt.begin(), basePt.end(), basePt.begin(),
                     std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0/lBar));
      air.densityFileIF = new DataFileIFReduced(inputFileNamePtr.c_str(), DataFileIFReduced::ASCII, IntVect(numGridPoints), RealVect(spacing), RealVect(basePt), 0, true);
      air.densityFileIF->GetAsciiData()->mult(1/nBar);
      air.densityFileIF->SetNoDataValue(air.densityFileIF->GetNoDataValue()/nBar);
      
    } else { // data from a function
      if (!pp.contains("densityProfilePieceNum")) {
        std::string profileName;
        int profileParamNum;
        vector<double> profileParam;
        pp.get("densityProfile", profileName);
        pp.get("densityProfileParamNum", profileParamNum);
        profileParam.resize(profileParamNum);
        pp.getarr("densityprofileParam", profileParam, 0, profileParamNum);
        if (profileName == "exp") {
          profileParam[0] = profileParam[0]/nBar;
          profileParam[1] = profileParam[1]/lBar;
          profileParam[2] = profileParam[2]/lBar;
        } else if (profileName == "stdAtm") {
          // need to do this because the function contains numbers with dimension
          profileParam[0] = profileParam[0];
          profileParam.push_back(lBar);
          profileParam.push_back(nBar);
        }
        
        air.bgdDensityProfile = new singleFunction(profileName, profileParam);
      } else {
        int numPiece;
        vector<double> lb;
        vector<string> fNames;
        vector<int>    paramNums;
        vector<vector<double>> paramVect;
        pp.get("densityProfilePieceNum", numPiece);
        pp.getarr("densityProfileLeftBound", lb, 0, numPiece);
        pp.getarr("densityProfileFuncNames", fNames, 0, numPiece);
        pp.getarr("densityProfileParamNums", paramNums, 0, numPiece);
        paramVect.resize(numPiece);
        for (int i = 0, startIdx = 0; i < numPiece; i++) {
          paramVect[i].resize(paramNums[i]);
          pp.getarr("densityProfileParamVect", paramVect[i], startIdx, paramNums[i]);
          startIdx += paramNums[i];
          if (fNames[i] == "exp") {
            lb[i] = lb[i]/lBar;
            paramVect[i][0] = paramVect[i][0]/nBar;
            paramVect[i][1] = paramVect[i][1]/lBar;
            paramVect[i][2] = paramVect[i][2]/lBar;
          }
        }
        air.bgdDensityProfile = new piecewiseFunction(numPiece, lb, fNames, paramVect);
      }
    }
  }
  
  
  std::string processName;
  int numPiece;
  vector<double> lb;
  vector<string> fNames;
  vector<int>    paramNums;
  vector<vector<double>> paramVect;
  
  pp = ParmParse("ionization");
  pp.get("name", processName);
  pp.get("numPieces", numPiece);
  lb.resize(numPiece);
  fNames.resize(numPiece);
  paramNums.resize(numPiece);
  pp.getarr("xlb", lb, 0, numPiece);
  pp.getarr("funcNames", fNames, 0, numPiece);
  pp.getarr("paramNums", paramNums, 0, numPiece);
  paramVect.resize(numPiece);
  for (int i = 0, startIdx = 0; i < numPiece; i++) {
    paramVect[i].resize(paramNums[i]);
    pp.getarr("A", paramVect[i], startIdx, paramNums[i]);
    startIdx += paramNums[i];
    lb[i] /= EBar/nBar;
    if (fNames[i] == "rcpExp") {
      paramVect[i][0] /= lBar*lBar;
      paramVect[i][1] /= EBar/nBar;
    } else if (fNames[i] == "const") {
      paramVect[i][0] /= lBar*lBar;
    } else {
      cout << "no coefficient nondimensionalization is done for process: " << processName << endl;
    }
  }
  
  air.EDpdentProcs[processName] = piecewiseFunction(numPiece, lb, fNames, paramVect);
  
  pp = ParmParse("attachment");
  pp.get("name", processName);
  pp.get("numPieces", numPiece);
  lb.resize(numPiece);
  fNames.resize(numPiece);
  paramNums.resize(numPiece);
  pp.getarr("xlb", lb, 0, numPiece);
  pp.getarr("funcNames", fNames, 0, numPiece);
  pp.getarr("paramNums", paramNums, 0, numPiece);
  paramVect.resize(numPiece);
  for (int i = 0, startIdx = 0; i < numPiece; i++) {
    paramVect[i].resize(paramNums[i]);
    pp.getarr("A", paramVect[i], startIdx, paramNums[i]);
    startIdx += paramNums[i];
    lb[i]    /= EBar/nBar;
    if (fNames[i] == "linear") {
      paramVect[i][0] /= lBar*lBar;
      paramVect[i][1] /= 1.0/(EBar*lBar);
    } else if (fNames[i] == "const") {
      paramVect[i][0] /= lBar*lBar;
    } else {
      cout << "no coefficient nondimensionalization is done for process: " << processName << endl;
    }
  }
  
  air.EDpdentProcs[processName] = piecewiseFunction(numPiece, lb, fNames, paramVect);
  
  pp = ParmParse("mobility");
  pp.get("name", processName);
  pp.get("numPieces", numPiece);
  lb.resize(numPiece);
  fNames.resize(numPiece);
  paramNums.resize(numPiece);
  pp.getarr("xlb", lb, 0, numPiece);
  pp.getarr("funcNames", fNames, 0, numPiece);
  pp.getarr("paramNums", paramNums, 0, numPiece);
  paramVect.resize(numPiece);
  for (int i = 0, startIdx = 0; i < numPiece; i++) {
    paramVect[i].resize(paramNums[i]);
    pp.getarr("A", paramVect[i], startIdx, paramNums[i]);
    startIdx += paramNums[i];
    lb[i] /= EBar/nBar;
    if (fNames[i] == "rcpLinear") {
      paramVect[i][0] /= muBar*nBar;
      paramVect[i][1] /= muBar*EBar;
    } else if (fNames[i] == "const") {
      paramVect[i][0] /= muBar*nBar;
    } else {
      cout << "no coefficient nondimensionalization is done for process: " << processName << endl;
    }
  }
  
  air.EDpdentProcs[processName] = piecewiseFunction(numPiece, lb, fNames, paramVect);
  
  a_gas = air;
  
  /*Real n;
   n = N;
   cout << "E=0: " << 1.0/lBar*n*a_gas.EDpdentProcs["ionization"].value(0.0) << "; " << "E=5: " << 1.0/lBar*n*a_gas.EDpdentProcs["ionization"].value(5/n) << "; " << "E=20: " << 1.0/lBar*n*a_gas.EDpdentProcs["ionization"].value(20/n) << endl;
   
   cout << "E=0: " << 1.0/lBar*n*a_gas.EDpdentProcs["attachment"].value(0.0) << "; " << "E=1.6: " << 1.0/lBar*n*a_gas.EDpdentProcs["attachment"].value(1.6/n) << "; " << "E=3.2: " << 1.0/lBar*n*a_gas.EDpdentProcs["attachment"].value(3.2/n) << "; " << "E=12.8: " << 1.0/lBar*n*a_gas.EDpdentProcs["attachment"].value(12.8/n) << endl;
   
   cout << "E=0.0: " << muBar/n*a_gas.EDpdentProcs["mobility"].value(0.0) << "; " << "E=0.002: " << muBar/n*a_gas.EDpdentProcs["mobility"].value(0.002/n) << "; " << "E=0.012: " << muBar/n*a_gas.EDpdentProcs["mobility"].value(0.012/n) << "; " << "E=0.2: " << muBar/n*a_gas.EDpdentProcs["mobility"].value(0.2/n) << "; " << "E=10: " << muBar/n*a_gas.EDpdentProcs["mobility"].value(2.0/n) << endl;
   
   cout << "n = " << n << endl;*/
  
}

void
getAdvectTestIBC(RefCountedPtr<AdvectTestIBC>& a_ibc)
{
  ParmParse pp("iniPlaCloud");
  
  MultiBlob blobs;

  if (pp.contains("randBlobNum")) {
    int blobNum;
    pp.get("randBlobNum", blobNum);
    if (blobNum > 0)
      generateRandomBlobs(blobs);
  }
  
  if (pp.contains("number")) {
    int blobNum;
    pp.get("number", blobNum);
    if (blobNum > 0)
      parseBlobsFromParmParse(blobs);
  }
  
  outputBlobs(blobs);
  
  pp = ParmParse("bgdPlasma");
  int numPiece;
  vector<double> lb;
  vector<string> fNames;
  vector<int>    paramNums;
  vector<vector<double>> paramVect;
  
  if (pp.contains("numPieces")) {
    
    pp.get("numPieces", numPiece);
    pp.getarr("xlb", lb, 0, numPiece);
    pp.getarr("funcNames", fNames, 0, numPiece);
    pp.getarr("paramNums", paramNums, 0, numPiece);
    paramVect.resize(numPiece);
    for (int i = 0, startIdx = 0; i < numPiece; i++) {
      paramVect[i].resize(paramNums[i]);
      pp.getarr("A", paramVect[i], startIdx, paramNums[i]);
      startIdx += paramNums[i];
      lb[i] /= lBar;
      if (fNames[i] == "Wait") {
        paramVect[i][0] /= nBar;
        paramVect[i][1] /= 1/lBar;
        paramVect[i][2] /= lBar;
        paramVect[i][3] /= 1/lBar;
        paramVect[i][4] /= lBar;
      } else if (fNames[i] == "tanhIP") {
        paramVect[i][0] /= nBar;
        paramVect[i][1] /= lBar;
        paramVect[i][2] /= lBar;
      } else if (fNames[i] == "const") {
        paramVect[i][0] /= nBar;
      } else {
        cout << "no parameter nondimensionalization is done for background plasma density profile!" << endl;
      }
    }
  } else {
    numPiece = 1;
    lb.push_back(-DBL_MAX);
    fNames.push_back("const");
    paramNums.push_back(1);
    paramVect.push_back(std::vector<double>{0.0});
  }
  
  piecewiseFunction bgdDensity = piecewiseFunction(numPiece, lb, fNames, paramVect);
  a_ibc = RefCountedPtr<AdvectTestIBC>(new AdvectTestIBC(blobs, bgdDensity));
}

void
getAMRLADFactory(RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
                 gas&                                          a_gas,
                 AdvectPhysics &                               a_advPhys)
{
  ParmParse pp;
  Real cfl = 0.8;
  pp.get("Ac",cfl);
  
  Real initialCFL = 0.1;
  pp.get("initial_Ac",initialCFL);
  
  Real domainLength;
  pp.get("domain_length",domainLength);
  domainLength = domainLength/normalization::scalingFactor/normalization::lBar;
  
  Real refineThresh = 0.3;
  pp.get ("refine_thresh",refineThresh);
  
  int tagBufferSize = 3;
  pp.get("tag_buffer_size",tagBufferSize);
  
  bool useLimiting = true;
  pp.get("use_limiting", useLimiting);
  
  Real nu;
  // pp.get("diffusion_coef", nu);
  // normalization DBar = lBar^2/tBar;
  // nu = nu / (normalization::lBar*normalization::lBar/normalization::tBar);
  nu = a_gas.m_elecDiffCoef;
  
  bool EPotBCVaryingFlag = false;
  ParmParse ppEPot("EPot");
  if (ppEPot.contains("bc_timeVarying"))
    ppEPot.get("bc_timeVarying", EPotBCVaryingFlag);
  
  if (!EPotBCVaryingFlag)
    a_fact = RefCountedPtr<AMRLevelAdvectDiffuseFactory>
    (new AMRLevelAdvectDiffuseFactory(a_advPhys, a_gas,
                                      ADParseBC, EPotParseBC, PIParseBC, cfl, domainLength,
                                      refineThresh, tagBufferSize,
                                      initialCFL, useLimiting, nu));
  else {
    RefCountedPtr<TimeDependentBCFunction> EPotBCFunc(new TimeDependentBCFunction());
    
    int numPiece;
    vector<double> lb;
    vector<string> fNames;
    vector<int>    paramNums;
    vector<vector<double>> uParamVect, udotParamVect;
    
    pp = ParmParse("bc_timeVarying");
    pp.get("numPieces", numPiece);
    lb.resize(numPiece);
    fNames.resize(numPiece);
    paramNums.resize(numPiece);
    pp.getarr("xlb", lb, 0, numPiece);
    pp.getarr("funcNames", fNames, 0, numPiece);
    pp.getarr("paramNums", paramNums, 0, numPiece);
    uParamVect.resize(numPiece);
    udotParamVect.resize(numPiece);
    for (int i = 0, startIdx = 0; i < numPiece; i++) {
      uParamVect[i].resize(paramNums[i]);
      pp.getarr("u", uParamVect[i], startIdx, paramNums[i]);
      udotParamVect[i].resize(paramNums[i]);
      pp.getarr("udot", udotParamVect[i], startIdx, paramNums[i]);
      startIdx += paramNums[i];
      lb[i] *= (1/scalingFactor) / tBar;
      if (fNames[i] == "linear") {
        uParamVect[i][0] /= phiBar;
        uParamVect[i][1] *= scalingFactor * (1/phiBar) * tBar; // V/s
        udotParamVect[i][0] *= scalingFactor / EBar;
        udotParamVect[i][1] *= scalingFactor * scalingFactor * (1/EBar) * tBar; // V/m/s
      } else if (fNames[i] == "const") {
        uParamVect[i][0] /= phiBar;
        udotParamVect[i][0] *= scalingFactor / EBar;
      } else {
        cout << "no coefficient nondimensionalization is done for time-varying BCs" << endl;
      }
    }
    
    if (procID() == uniqueProc(SerialTask::compute)) {
      pout() << "time-varying BCs parameters" << endl;
      pout() << "functions = ";
      for (const auto& fname : fNames) {
        pout() << fname << " ";
      }
      pout() << endl;
      pout() << "lb = ";
      for (const auto& lbvalue : lb) {
        pout() << lbvalue*tBar << " ";
      }
      pout() << endl;
      pout() << "uParam = " << uParamVect[0][0]*phiBar << " " << uParamVect[0][1]*phiBar*(1/tBar) << "; " << uParamVect[1][0]*phiBar << endl;
      pout() << "udotParam = " << udotParamVect[0][0]*EBar << " " << udotParamVect[0][1]*EBar*(1/tBar) << "; " << udotParamVect[1][0]*EBar << endl;
    }
    // Now create two piecewise functions: one for , one for udot
    piecewiseFunction u(numPiece, lb, fNames, uParamVect);
    piecewiseFunction udot(numPiece, lb, fNames, udotParamVect);
    
    // Assign to EPotBCFunc
    EPotBCFunc->m_valueFunc->m_timeFunctions.clear();
    EPotBCFunc->m_valueFunc->m_timeFunctions.push_back(u);
    EPotBCFunc->m_valueFunc->m_timeFunctions.push_back(udot);
    
    EPotBCFunc->setValueFunction(linearTimeLinearZ);
    
    BCHolder EPotTimeVaryBC(EPotBCFunc);
    a_fact = RefCountedPtr<AMRLevelAdvectDiffuseFactory>(
                                                         new AMRLevelAdvectDiffuseFactory(a_advPhys, a_gas,
                                                                                          ADParseBC, EPotTimeVaryBC, PIParseBC, cfl, domainLength,
                                                                                          refineThresh, tagBufferSize,
                                                                                          initialCFL, useLimiting, nu));
  }
}

void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat)
{
  ParmParse pp;
  int max_level = 0;
  pp.get("max_level",max_level);
  
  int num_read_levels = Max(max_level,1);
  std::vector<int> regrid_intervals; // (num_read_levels,1);
  pp.getarr("regrid_interval",regrid_intervals,0,num_read_levels);
  
  int block_factor = 1;
  pp.get("block_factor",block_factor);
  
  int max_grid_size = 32;
  pp.get("max_grid_size",max_grid_size);
  
  Real fill_ratio = 0.75;
  pp.get("fill_ratio",fill_ratio);
  
  int checkpoint_interval = 0;
  pp.get("checkpoint_interval",checkpoint_interval);
  
  int plot_interval = 0;
  pp.get("plot_interval",plot_interval);
  
  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  pp.query("fixed_dt",fixedDt);
  
  Real max_dt_growth = 1.1;
  pp.get("max_dt_growth",max_dt_growth);
  
  Real dt_tolerance_factor = 1.1;
  pp.get("dt_tolerance_factor",dt_tolerance_factor);
  a_amr.define(max_level, a_refRat,
               a_prob_domain,&(*a_fact));
  
  if (fixedDt > 0)
  {
    a_amr.fixedDt(fixedDt/normalization::tBar);
  }
  
  // set grid generation parameters
  a_amr.maxGridSize(max_grid_size);
  a_amr.blockFactor(block_factor);
  a_amr.fillRatio(fill_ratio);
  
  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize;
  pp.get("grid_buffer_size",gridBufferSize);
  a_amr.gridBufferSize(gridBufferSize);
  
  // set output parameters
  a_amr.checkpointInterval(checkpoint_interval);
  a_amr.plotInterval(plot_interval);
  a_amr.regridIntervals(regrid_intervals);
  a_amr.maxDtGrow(max_dt_growth);
  a_amr.dtToleranceFactor(dt_tolerance_factor);
  
  if (pp.contains("plot_period")) {
    double plot_period = 0.0;
    pp.get("plot_period", plot_period);
    a_amr.plotPeriod(plot_period/normalization::scalingFactor/normalization::tBar);
  }
  
  if (pp.contains("use_subcycling"))
  {
    bool useSubcycling;
    pp.get("use_subcycling", useSubcycling);
    if (!useSubcycling)
    {
      pout() << "SUBCYCLING IN TIME TURNED OFF!!!"  << endl;
    }
    a_amr.useSubcyclingInTime(useSubcycling);
  }
  if (pp.contains("plot_prefix"))
  {
    std::string prefix;
    pp.get("plot_prefix",prefix);
    a_amr.plotPrefix(prefix);
    //    pout() << procID() << " " << uniqueProc(SerialTask::compute) << endl;
    if (procID() == uniqueProc(SerialTask::compute)) {
      std::size_t found = prefix.find_last_of("/");
      if (found!=std::string::npos) {
        struct stat info;
        std::string dirNameString = prefix.substr(0, found);
        char* dirName = (char*) (dirNameString.c_str());
        int status = stat(dirName, &info);
        if(status < 0) {
          pout() << "Output folder doesn't exist, and one is created" << endl;
          mkdir(dirName, 0777);
        } else if (info.st_mode & S_IFDIR)
          pout() << "Use existing output folder: " << dirNameString << endl;
        else
          pout() << "Invalid output folder name!" << endl;
      }
    }
  }
  
  if (pp.contains("chk_prefix"))
  {
    std::string prefix;
    pp.get("chk_prefix",prefix);
    a_amr.checkpointPrefix(prefix);
  }
  
  int verbosity;
  pp.get("verbosity",verbosity);
  CH_assert(verbosity >= 0);
  
  a_amr.verbosity(verbosity);
}
void
setupAMRForAMRRun(AMR& a_amr)
{
  ParmParse pp;
  
  if (!pp.contains("restart_file"))
  {
    // initialize from scratch for AMR run
    // initialize hierarchy of levels
    a_amr.setupForNewAMRRun();
  }
  else
  {
    std::string restart_file;
    pp.get("restart_file",restart_file);
    pout() << " restarting from file " << restart_file << endl;
    
#ifdef CH_USE_HDF5
    HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    a_amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
  }
  
}

void averageCellToEdge (const DisjointBoxLayout a_grids, LevelData<FArrayBox>& a_U, LevelData<FluxBox>& a_UEdge) {
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit) {
    //    Box b = a_U[dit].box();
    Box b = a_grids[dit];
    FArrayBox& UFab = a_U[dit];
    
    for (int dir = 0; dir < SpaceDim; ++dir) {
      const Box bCenter = b & grow(a_grids.physDomain(), -BASISV(dir));
      // const Box bCenter1 = b & a_grids.physDomain();
      const Box bLo     = b & adjCellLo(bCenter,dir);
      const Box bHi     = b & adjCellHi(bCenter,dir);
      //      const Box bLoGhost     = b & adjCellLo(a_grids.physDomain(),dir);
      //      const Box bHiGhost     = b & adjCellHi(a_grids.physDomain(),dir);
      const Box bLoGhost = adjCellBox(b, dir, Side::Lo, 1);
      const Box bHiGhost = adjCellBox(b, dir, Side::Hi, 1);
      // copy the boundary cells of the orignal box
      if (!bLo.isEmpty()) UFab.copy(UFab, bLo, 0, bLoGhost, 0, a_U.nComp());
      if (!bHi.isEmpty()) UFab.copy(UFab, bHi, 0, bHiGhost, 0, a_U.nComp());
    }
  }
  
  
  //  printDiagnosticInfo (100, 1, a_grids,a_U, "advVelCC", "averageCellToEdge");
  //  printDiagnosticInfo (100, 1, a_grids,a_UEdge, "advVel", "averageCellToEdge");
  
  CellToEdge(a_U, a_UEdge);
  
  //  printDiagnosticInfo (100, 1, a_grids,a_UEdge, "advVel", "averageCellToEdge");
  
  //  double sum = 0;
  //  double FluxMin;
  //  double FluxMax;
  //  IntVect minLoc(-100, -100, -100), maxLoc(-100, -100, -100);
  //
  //  int comp = 0;
  //  for (int dir = 0; dir < SpaceDim; ++dir) {
  //    sum = 0;
  //    FluxMin = numeric_limits<double>::max();
  //    FluxMax = numeric_limits<double>::lowest();
  //    minLoc = IntVect(-100, -100, -100);
  //    maxLoc = IntVect(-100, -100, -100);
  //    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++ dit) {
  //      const FArrayBox& fab = a_UEdge[dit()].getFlux(dir);
  ////      const Box& b = a_grids[dit()].surroundingNodes(dir);
  //      const Box& b = fab.box(); //this includes ghost cells
  //      sum += fab.sum(b, comp);
  //      if(fab.min(b, comp) < FluxMin) {
  //        FluxMin = fab.min(b, comp);
  //        minLoc = fab.minIndex(b, comp);
  //      }
  //      if(fab.max(b, comp) > FluxMax) {
  //        FluxMax = fab.max(b, comp);
  //        maxLoc = fab.maxIndex(b, comp);
  //      }
  //    }
  //  }
  //
  //  if (FluxMax > 1e100)
  //    pout() << "Error!";
  
}

void unnormalize (LevelData<FArrayBox>& a_U, int a_comp, int a_numcomp, Real norm) {
  DisjointBoxLayout grids = a_U.disjointBoxLayout();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit) {
    FArrayBox& UFab = a_U[dit];
    UFab.mult(norm, a_comp, a_numcomp);
  }
}

/*******/
void printDiagnosticInfo (int a_level, double a_dx, DisjointBoxLayout a_grids, const LevelData<FArrayBox>& U, const std::string& variableName, const std::string& funcInfo) {
  
  pout() << left << setw(50) << funcInfo << " L" << a_level << endl;
  
  double sum = 0;
  double UMin = numeric_limits<double>::max();
  double UMax = numeric_limits<double>::lowest();
  IntVect minLoc(-100, -100, -100), maxLoc(-100, -100, -100);
  for (int comp = 0; comp < U.nComp(); comp++) {
    sum = 0;
    UMin = numeric_limits<double>::max();
    UMax = numeric_limits<double>::lowest();
    minLoc = IntVect(-100, -100, -100);
    maxLoc = IntVect(-100, -100, -100);
    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++ dit) {
      //      const Box& b = a_grids[dit()];
      const Box& b = U[dit()].box();
      sum += U[dit()].sum(b, comp);
      if(U[dit()].min(b, comp) < UMin) {
        UMin = min(U[dit()].min(b, comp), UMin);
        minLoc = U[dit()].minIndex(b, comp);
      }
      if(U[dit()].max(b, comp) > UMax) {
        UMax = max(U[dit()].max(b, comp), UMax);
        maxLoc = U[dit()].maxIndex(b, comp);
      }
    }
    // this is to get the total particles if U represents density
    sum *= pow(a_dx, SpaceDim);
    
    std::string smin, smax;
    std::stringstream stros;
    stros.str(std::string());
    stros.clear();
    stros << minLoc;
    smin = stros.str();
    smin.append(20-smin.length(), ' ');
    
    stros.str(std::string());
    stros.clear();
    stros << maxLoc;
    smax = stros.str();
    smax.append(20-smax.length(), ' ');
    
    pout() << left << setw(10) << variableName << " com " << comp << " min " << "@" << left << smin << left << setw(16) << UMin << " max " << "@" << left << smax << left << setw(16) << UMax << "  sum " << left << setw(16) << sum << endl;
  }
}

/*******/
void printDiagnosticInfo (int a_level, double a_dx, DisjointBoxLayout a_grids, const LevelData<FluxBox>& a_levelFlux, const std::string& variableName, const std::string& funcInfo) {
  
  double sum = 0;
  double FluxMin;
  double FluxMax;
  IntVect minLoc(-100, -100, -100), maxLoc(-100, -100, -100);
  
  pout() << left << setw(50) << funcInfo << " L" << a_level << endl;
  
  int comp = 0;
  for (int dir = 0; dir < SpaceDim; ++dir) {
    sum = 0;
    FluxMin = numeric_limits<double>::max();
    FluxMax = numeric_limits<double>::lowest();
    minLoc = IntVect(-100, -100, -100);
    maxLoc = IntVect(-100, -100, -100);
    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++ dit) {
      const FArrayBox& fab = a_levelFlux[dit()].getFlux(dir);
      //Box b = a_grids[dit()].surroundingNodes(dir);
            const Box& b = fab.box(); //this includes ghost cells
      //      cout << b.type() << " " << b << endl;
      //      cout << fab.box().type() << " " << fab.box() << endl;
      sum += fab.sum(b, comp);
      if(fab.min(b, comp) < FluxMin) {
        FluxMin = fab.min(b, comp);
        minLoc = fab.minIndex(b, comp);
      }
      if(fab.max(b, comp) > FluxMax) {
        FluxMax = fab.max(b, comp);
        maxLoc = fab.maxIndex(b, comp);
      }
      
      //      for (BoxIterator bit(b); bit.ok(); ++bit) {
      //        const IntVect& iv = bit();
      //        if(FluxMin > fab(iv, 0)) {
      //          FluxMin = fab(iv, 0);
      //          minLoc = bit();
      //        }
      //
      //        if(FluxMax < fab(iv, 0)) {
      //          FluxMax = fab(iv, 0);
      //          maxLoc = bit();
      //        }
      //      }
    }
    
    // this is to get the total particles if U represents density
    sum *= pow(a_dx, SpaceDim);
    
    std::stringstream stros;
    std::string smin, smax;
    
    stros.str(std::string());
    stros << minLoc;
    smin = stros.str();
    smin.append(20-smin.length(), ' ');
    
    stros.str(std::string());
    stros << maxLoc;
    smax = stros.str();
    smax.append(20-smax.length(), ' ');
    
    pout() << left << setw(10) << variableName << " dir " << dir << " min " << "@" << left << smin << left << setw(16) << FluxMin << " max " << "@" << left << smax << left << setw(16) << FluxMax << "   sum " << left << setw(16) << sum << endl;
  }
}

/*******/
void outputDataForCheck (int a_level, DisjointBoxLayout a_grids, const LevelData<FluxBox>& a_levelFlux) {
  pout() << endl << "------------ level " << a_level << " ---------------";
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit) {
    pout() << endl << "------------ level " << a_level << " -- box: ";
    dit.i().output();
    for (int dir = 0; dir < SpaceDim; ++dir) {
      pout() << endl << "------------ level " << a_level << " -- dir " << dir << " ---------------" << endl;
      const FArrayBox& fab = a_levelFlux[dit()].getFlux(dir);
      //      FArrayBox& fab = a_levelData[dit()];
      const Box& fabBox = fab.box();
      for (BoxIterator bit(fabBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        pout() << iv << ":" << fab(iv, 0) << "  ";
      }
    }
  }
  pout() << endl;
}

/*******/
void outputDataForCheck (int a_level, DisjointBoxLayout a_grids, const LevelData<FArrayBox>& a_U) {
  pout() << endl << "------------------ level " << a_level << " -------------------";
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit) {
    pout() << endl << "------------ level " << a_level << " -- box: ";
    dit.i().output();
    for (int iComp = 0; iComp < a_U.nComp(); ++iComp) {
      pout() << endl << "------------ level " << a_level << " -- comp " << iComp << " ---------------" << endl;
      const FArrayBox& fab = a_U[dit()];
      //      FArrayBox& fab = a_levelData[dit()];
      const Box& fabBox = fab.box();
      for (BoxIterator bit(fabBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        pout() << iv << ":" << fab(iv, iComp) << "  ";
      }
    }
  }
  pout() << endl;
}

#include "NamespaceFooter.H"
