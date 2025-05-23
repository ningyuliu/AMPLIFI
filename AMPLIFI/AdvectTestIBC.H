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

#ifndef _ADVECTTESTIBC_H_
#define _ADVECTTESTIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysIBC.H"
#include "SolidBCF_F.H"
#include "parameterizedFunction.hpp"
#include "blob.hpp"

#include "NamespaceHeader.H"

/// Example IBC for advection-diffusion problem
class AdvectTestIBC:public PhysIBC
{
public:

  ///
  AdvectTestIBC(const MultiBlob&             a_blob,
                const piecewiseFunction&     a_bgdDensity)
  {
    m_blob       = a_blob;
    m_bgdDensity = a_bgdDensity;
  }
  
  /// Destructor
  ~AdvectTestIBC()
  {
  };

  /// Factory method - this object is its own factory
  PhysIBC *new_physIBC();

  /// Set boundary fluxes
  void primBC(FArrayBox&            a_F,
              const FArrayBox&      a_W,
              const FArrayBox&      a_Wextrap,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  ///one sided is OK by me
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// does not apply here
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time)
  {
  }

  /// Set up initial conditions
  void initialize(LevelData<FArrayBox>& a_U);

protected:

  MultiBlob           m_blob;
  piecewiseFunction   m_bgdDensity;

private:
  //weak construction is bad
  AdvectTestIBC()
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  void operator=(const AdvectTestIBC&);
  AdvectTestIBC(const AdvectTestIBC&);

};

#include "NamespaceFooter.H"

#endif

