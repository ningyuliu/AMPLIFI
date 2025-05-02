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
#include "hyperbolicBC.hpp"

#include "LoHiCenter.H"
#include "NamespaceHeader.H"

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used).
PhysIBC* AdvectTestIBC::new_physIBC()
{
  AdvectTestIBC* retval = new AdvectTestIBC(m_number, m_center, m_a, m_mag, m_blob, m_bgdDensity);
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
      
      RealVect point(iv);
      RealVect ccOffset = 0.5*m_dx*RealVect::Unit;
      point *= m_dx;
      point += ccOffset;
      vector<double> pvec(SpaceDim);
      a_U[dit()](iv, 0) = m_blob.value(pvec);
      for(int i = 0; i != pvec.size(); i++)
        pvec[i] = point[i];
      a_U[dit()](iv, 0) += m_bgdDensity.value(pvec);
    }    
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
//          FORT_SOLIDBCF(CHF_FRA(a_WGdnv),
//                        CHF_CONST_FRA(a_Wextrap),
//                        CHF_CONST_FRA(a_W),
//                        CHF_CONST_INT(lohisign),
//                        CHF_CONST_INT(a_dir),
//                        CHF_BOX(boundaryBox));
          
          SolidBC(a_WGdnv, a_Wextrap, a_W, a_dir, lohisign, boundaryBox);
        }
    }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void AdvectTestIBC::setBdrySlopes(FArrayBox&       a_dW,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Real&      a_time)
{
  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    Box loBox,hiBox,centerBox,domain;
    int hasLo,hasHi;
    Box slopeBox = a_dW.box();
    slopeBox.grow(a_dir,1);

    // Generate the domain boundary boxes, loBox and hiBox, if there are
    // domain boundarys there
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
               slopeBox,m_domain,a_dir);

    // Set the boundary slopes if necessary
    if ((hasLo != 0) || (hasHi != 0))
    {
//      FORT_SLOPEBCSF(CHF_FRA(a_dW),
//                     CHF_CONST_FRA(a_W),
//                     CHF_CONST_INT(a_dir),
//                     CHF_BOX(loBox),
//                     CHF_CONST_INT(hasLo),
//                     CHF_BOX(hiBox),
//                     CHF_CONST_INT(hasHi));
      
      SlopeBC(a_dW, a_W, a_dir, loBox, hasLo, hiBox, hasHi);
    }
  }
}

#include "NamespaceFooter.H"
