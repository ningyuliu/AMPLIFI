//
//  hyperbolicBC.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 1/2/25.
//

#ifndef hyperbolicBC_hpp
#define hyperbolicBC_hpp

#include <stdio.h>

#include "FArrayBox.H"

namespace hyperbolicBC {
  
  const double smallDensity = 1e-10;
  
}

void SolidBC(FArrayBox&            a_WGdnv,
             const FArrayBox&      a_Wextrap,
             const FArrayBox&      a_W,
             const int&            a_dir,
             const int&            isign,
             const Box&            boundaryBox);

void SlopeBC(FArrayBox&            a_dW,
             const FArrayBox&      a_W,
             const int&            a_dir,
             const Box&            loBox,
             const int&            hasLo,
             const Box&            hiBox,
             const int&            hasHi);

#endif /* hyperbolicBC_hpp */
