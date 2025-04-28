//
//  additionalBCFunc.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#include "additionalBCFunc.hpp"

#include "BCFunc.H"
#include "RealVect.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

//The following function is from Chombo BCFunc.cpp
static void getDomainFacePosition(RealVect&             a_retval,
                                  const IntVect&        a_validIV,
                                  const Real&           a_dx,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side)
{
  Real* dataPtr = a_retval.dataPtr();

  D_TERM( dataPtr[0] = a_dx*(a_validIV[0] + 0.5);,\
          dataPtr[1] = a_dx*(a_validIV[1] + 0.5);,\
          dataPtr[2] = a_dx*(a_validIV[2] + 0.5);)

  int isign = sign(a_side);
  dataPtr[a_dir] += 0.5*Real(isign)*a_dx;
}

// Implement Robin BC: (alpha + d/dx) phi = gamma
///
/**
 Robin  boundary conditions for a side.
 */
void RobinBC(FArrayBox& a_state, const Box& a_valid, Real a_dx, bool a_homogeneous, const BCValueHolder& a_valueA, int a_dir, Side::LoHiSide a_side) {
  
  BCValueHolder& a_value = (BCValueHolder&)a_valueA;
  int isign = sign(a_side);
  RealVect facePos;
  
  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();
  
  Real* value = new Real[2];
  for (BoxIterator bit(toRegion); bit.ok(); ++bit) {
    const IntVect& ivTo = bit();
    IntVect ivClose = ivTo -   isign*BASISV(a_dir);
    Real nearVal = a_state(ivClose);
    Real alpha, gamma;

    getDomainFacePosition(facePos, ivClose, a_dx, a_dir, a_side);
    a_value(facePos.dataPtr(), &a_dir, &a_side, value);
    
    alpha = value[0];
    gamma = value[1];
    a_state(ivTo) = (1 - Real(isign)*a_dx*alpha) * nearVal + Real(isign)*a_dx*gamma;
  }
  
  delete[] value;
}

#include "NamespaceFooter.H"
