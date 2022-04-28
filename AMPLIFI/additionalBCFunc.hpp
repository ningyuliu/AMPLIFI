//
//  additionalBCFunc.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef additionalBCFunc_hpp
#define additionalBCFunc_hpp

// implementation of general BCs that are unavailable in Chombo
#include <stdio.h>
#include "BCFunc.H"

void RobinBC (FArrayBox& a_state, const Box& a_valid, Real a_dx, bool a_homogeneous, const BCValueHolder& a_value, int a_dir, Side::LoHiSide a_side);

#endif /* additionalBCFunc_hpp */
