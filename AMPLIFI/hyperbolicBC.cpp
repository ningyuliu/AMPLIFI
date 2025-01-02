//
//  hyperbolicBC.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 1/2/25.
//

#include "hyperbolicBC.hpp"

void SolidBC(FArrayBox&            a_WGdnv,
             const FArrayBox&      a_Wextrap,
             const FArrayBox&      a_W,
             const int&            a_dir,
             const int&            isign,
             const Box&            boundaryBox) {
    
  for (BoxIterator bit(boundaryBox); bit.ok(); ++bit) {
    const IntVect& ivTo = bit();
    a_WGdnv(ivTo) = max(a_Wextrap(ivTo), hyperbolicBC::smallDensity);
  }
}

void SlopeBC(FArrayBox&            a_dW,
             const FArrayBox&      a_W,
             const int&            a_dir,
             const Box&            loBox,
             const int&            hasLo,
             const Box&            hiBox,
             const int&            hasHi) {
  
  if (hasLo)
    for (BoxIterator bit(loBox); bit.ok(); ++bit) {
      const IntVect& ivTo = bit();
      a_dW(ivTo) = 0.0;
    }
  
  if (hasHi)
    for (BoxIterator bit(hiBox); bit.ok(); ++bit) {
      const IntVect& ivTo = bit();
      a_dW(ivTo) = 0.0;
    }
}
