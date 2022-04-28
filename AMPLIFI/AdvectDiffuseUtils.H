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

#ifndef _ADVECTDIFFUSEUTILS_H_
#define _ADVECTDIFFUSEUTILS_H_

#include <cmath>

#include "AMR.H"
#include "AMRLevelAdvectDiffuseFactory.H"
#include "AdvectTestIBC.H"
#include "AdvectionFunctions.H"
#include "memtrack.H"
#include "parstream.H"

#include "NamespaceHeader.H"

///
/**
 */
extern void
makeFinestDomain(ProblemDomain& a_domain,
                 Real&          a_dx);

///
/**
**/
extern Real
dtgNorm(const Vector< LevelData<FArrayBox>* >& a_src,
        const Vector< DisjointBoxLayout >&     a_grids,
        const Vector<int>&                     a_refRatio,
        const ProblemDomain&                   a_coarseDom,
        const int& a_comp,
        const int& a_pval);

///
/**
**/
extern void
compareError(Vector<Real>&                            a_orders,
             const Vector< LevelData<FArrayBox>* >&   a_errorFine,
             const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
             const Vector< DisjointBoxLayout >&       a_gridsFine,
             const Vector< DisjointBoxLayout >&       a_gridsCoar,
             const Vector<int>&                       a_refRat,
             const ProblemDomain&                     a_coarseDom,
             int a_testverbosity);


///
/**
 */
extern void
getErrorFromCoarseAndFine(Vector< LevelData<FArrayBox>* >&           a_errorCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnCoar,
                          const Vector< DisjointBoxLayout >&         a_gridsCoar,
                          const ProblemDomain&                       a_level0DomainCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnFine,
                          const Vector< DisjointBoxLayout >&         a_gridsFine,
                          const ProblemDomain&                       a_level0DomainFine,
                          const Vector<int>&                         a_refRat);
///
/**
**/
extern void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar);
///
/*****/
extern void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain);

///
/**
 */
extern void
getProblemDomain(ProblemDomain& a_domain);

///
/**
 */
extern void
getGlobalVariables();

///
/**
 */
extern void
getAmbientGas(gas& a_gas);

///
/**
 */
extern void
getAdvectTestIBC(RefCountedPtr<AdvectTestIBC>&  a_advPhys);

///
/**
 */
extern void
getAMRLADFactory(RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
                 gas&                                          a_gas,
                 AdvectPhysics &                               a_advPhys);


///
/**
 */
extern void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat);

///
/**
 */
extern void
setupAMRForAMRRun(AMR& a_amr);

/// Average cell data to edges
/** Ghost cells are simply filled with a copy of the cells next to the boundary
 */
extern void
averageCellToEdge (const DisjointBoxLayout a_grids, LevelData<FArrayBox>& a_U, LevelData<FluxBox>& a_UEdge);

/// Unnormalize the variables for outputting
/** a_U = a_U * norm
 */
extern void
unnormalize (LevelData<FArrayBox>& a_U, int a_comp, int a_numcomp, Real norm);

/// print the min, max and sum of the part of U on each processor. Future work: print these values from the
/// entire U.
/**
 */
extern void
printDiagnosticInfo (int a_level, double a_dx, DisjointBoxLayout a_grids, const LevelData<FArrayBox>& U, const std::string& variableName, const std::string& funcInfo);

/// print the min, max and sum of the part of U on each processor. Future work: print these values from the
/// entire U.
/**
 */
extern void
printDiagnosticInfo (int a_level, double a_dx, DisjointBoxLayout a_grids, const LevelData<FluxBox>& a_levelFlux, const std::string& variableName, const std::string& funcInfo);

/// output data for checking
/**
 */
extern void
outputDataForCheck (int a_Level, DisjointBoxLayout a_grids, const LevelData<FluxBox>& a_levelFlux);

/// output data for checking
/**
 */
void outputDataForCheck (int a_level, DisjointBoxLayout a_grids, const LevelData<FArrayBox>& a_U);


/****/


#include "NamespaceFooter.H"

#endif
