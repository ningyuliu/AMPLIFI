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

#include "LoHiSide.H"

#include "AdvectTestIBC.H"

#include "LoHiCenter.H"
#include "NamespaceHeader.H"

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used).
PhysIBC* AdvectTestIBC::new_physIBC()
{
  AdvectTestIBC* retval = new AdvectTestIBC(m_number, m_center, m_a, m_mag, m_bgdDensity);
  return static_cast<PhysIBC*>(retval);
}

// Set up initial conditions
void AdvectTestIBC::initialize(LevelData<FArrayBox>& a_U)
{
  DisjointBoxLayout grids = a_U.disjointBoxLayout();
  LevelData<FArrayBox> tmp;
  tmp.define(a_U);
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
  {
    a_U[dit()].setVal(0.0, a_U[dit()].box(), 0);
    tmp[dit()].setVal(0.0, tmp[dit()].box(), 0);
    
//    for (int i = 0; i < m_number; i++) {
//      FORT_ADVECTINITF(CHF_FRA1(tmp[dit()],0),
//                       CHF_CONST_REALVECT(m_center[i]),
//                       CHF_CONST_REAL(m_a[i]),
//                       CHF_CONST_REAL(m_dx),
//                       CHF_BOX(grid));
//      tmp[dit()].mult(m_mag[i]);
//      a_U[dit()].plus(tmp[dit()]);
//    }
    
    for (BoxIterator bit(grids.get(dit)); bit.ok(); ++bit) {
      const IntVect& iv = bit();
      for (int i = 0; i < m_number; i++) {
        double val = m_mag[i];
        for (int dir = 0; dir < SpaceDim; dir++) {
          double tmp = 0;
          tmp = (iv[dir]+0.5)*m_dx - m_center[i][dir];
          tmp = tmp*tmp/m_a[i][dir]/m_a[i][dir];
          val *= exp(-tmp);
        }
        a_U[dit()](iv, 0) += val;
      }
    }
    a_U[dit()] += m_bgdDensity;
    
  }
}

// Set boundary fluxes
void AdvectTestIBC::primBC(FArrayBox&            a_WGdnv,
                           const FArrayBox&      a_Wextrap,
                           const FArrayBox&      a_W,
                           const int&            a_dir,
                           const Side::LoHiSide& a_side,
                           const Real&           a_time)
{
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
    {
      int lohisign;
      Box tmp = a_WGdnv.box();

      // Determine which side and thus shifting directions
      lohisign = sign(a_side);
      tmp.shiftHalf(a_dir,lohisign);

      // Is there a domain boundary next to this grid
      if (!m_domain.contains(tmp))
        {
          tmp &= m_domain;

          Box boundaryBox;

          // Find the strip of cells next to the domain boundary
          if (a_side == Side::Lo)
            {
              boundaryBox = bdryLo(tmp,a_dir);
            }
          else
            {
              boundaryBox = bdryHi(tmp,a_dir);
            }

          // Set the boundary fluxes
          FORT_SOLIDBCF(CHF_FRA(a_WGdnv),
                        CHF_CONST_FRA(a_Wextrap),
                        CHF_CONST_FRA(a_W),
                        CHF_CONST_INT(lohisign),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(boundaryBox));
        }
    }
}

#include "NamespaceFooter.H"
