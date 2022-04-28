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

#ifndef _AMRLEVELADVECTDIFFUSEFACTORY_H_
#define _AMRLEVELADVECTDIFFUSEFACTORY_H_

#include "AMRLevelFactory.H"

#include "AdvectionFunctions.H"
#include "AMRLevelAdvectDiffuse.H"

#include "NamespaceHeader.H"

///
/**
 */
class AMRLevelAdvectDiffuseFactory: public AMRLevelFactory
{
public:

  /// Constructor
  /**
   */
  AMRLevelAdvectDiffuseFactory(/// advection physics class
                               const AdvectPhysics&        a_gphys,
                               gas&                        a_gas,
                               /// boundary condition class for diffusion solve
                               BCHolder                    a_bcFunc,
                               /// boundary condition class for field solve
                               BCHolder                    a_EPotbcFunc,
                               /// boundary condition class for photoionization solve
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


  /// Destructor
  virtual ~AMRLevelAdvectDiffuseFactory()
  {
  }


  ///
  virtual AMRLevel* new_amrlevel() const;

protected:

  Real                  m_cfl;
  Real                  m_domainLength;
  Real                  m_refineThresh;
  int                   m_tagBufferSize;
  Real                  m_initialDtMultiplier;
  bool                  m_useLimiting;
  Real                  m_nu;
  RefCountedPtr<AdvectPhysics>     m_advPhys;
  gas                   m_gas;
  BCHolder              m_bcFunc; //bc for diffusion solve
  BCHolder              m_EPotbcFunc; //bc for field solve
  BCHolder              m_phtznbcFunc; //bc for field solve

private:

  /// Weak construction is disallowed
  AMRLevelAdvectDiffuseFactory()
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelAdvectDiffuseFactory&);
  AMRLevelAdvectDiffuseFactory(const AMRLevelAdvectDiffuseFactory&);
};

#include "NamespaceFooter.H"

#endif
