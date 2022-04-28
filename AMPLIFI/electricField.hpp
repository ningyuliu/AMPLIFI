//
//  electricField.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef electricField_hpp
#define electricField_hpp

// This is for representing electric field and related boundary conditions

#include "RealVect.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

// define a BC with constant value
extern void PoissonBCParseValue(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);

// const BC at zlow, zhigh, but linearly varies with z along the other
// boundaries, where bcValLo/Hi stores the -rate/field
extern void PoissonBCLinearAlongZ(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);

extern void EPotParseBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous);

class field {
  
public:
  field() {}
  void define(DisjointBoxLayout a_grids, IntVect ivGhost1, IntVect ivGhost2);
  void define(const field& a_field);
  void copyTo(field& des);
  virtual ~field() {};
  
  LevelData<FArrayBox> m_E;             // E vector at cell center
  LevelData<FArrayBox> m_Emag;          // E magnitude at cell center
  LevelData<FluxBox>   m_EEdge;         // Edge field
  
protected:
  // verbosity level
  static int s_verbosity;
  
};


#include "NamespaceFooter.H"

#endif /* electricField_hpp */
