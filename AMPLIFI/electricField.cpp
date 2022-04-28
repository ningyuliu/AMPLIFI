//
//  electricField.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#include "ParmParse.H"
#include "RealVect.H"
#include "IntVect.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LayoutIterator.H"
#include "BCFunc.H"

#include "electricField.hpp"

#include "NamespaceHeader.H"

int field::s_verbosity;

void PoissonBCParseValue(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values) {

  ParmParse pp("EPot");
  std::vector<Real> bcVal;

  if(*side == Side::Lo)
    pp.getarr("bc_loValue", bcVal, 0, SpaceDim);
  else
    pp.getarr("bc_hiValue", bcVal, 0, SpaceDim);

  a_values[0] = bcVal[*dir];
}

// const. at zlow, zhigh, but linearly varies with z along the other boundaries, where bcValLo/Hi stores the -rate/field
void PoissonBCLinearAlongZ(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values) {
  
  ParmParse pp("EPot");
  std::vector<Real> bcValLo, bcValHi;
  
  pp.getarr("bc_loValue", bcValLo, 0, SpaceDim);
  pp.getarr("bc_hiValue", bcValHi, 0, SpaceDim);
  
  if (*dir == SpaceDim-1)
    if(*side == Side::Lo)
      a_values[0] = bcValLo[SpaceDim-1];
    else
      a_values[0] = bcValHi[SpaceDim-1];
  else
    if(*side == Side::Lo)
      a_values[0] = bcValLo[SpaceDim-1] - pos[SpaceDim-1]*bcValLo[*dir];
    else
      a_values[0] = bcValLo[SpaceDim-1] - pos[SpaceDim-1]*bcValHi[*dir];

}

void EPotParseBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous) {
  
  if (!a_domain.domainBox().contains(a_state.box())) {
    
    Box valid = a_valid;
    int a_order = 1;
    
    for (int i=0; i<CH_SPACEDIM; ++i)
      
      // don't do anything if periodic
      if (!a_domain.isPeriodic(i)) {
        ParmParse pp("EPot");
        std::vector<int>  bcLo = std::vector<int>();
        std::vector<int>  bcHi = std::vector<int>();
        pp.getarr("bc_lo", bcLo, 0, SpaceDim);
        pp.getarr("bc_hi", bcHi, 0, SpaceDim);
        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        
        if (!a_domain.domainBox().contains(ghostBoxLo))
          
          switch (bcLo[i]) {
            case 0:
              DiriBC(a_state, valid, a_dx, a_homogeneous, PoissonBCParseValue, i, Side::Lo);
              break;
              
            case 1:
              NeumBC(a_state, valid, a_dx, a_homogeneous, PoissonBCParseValue, i, Side::Lo);
              break;
              
            case 2:
              ExtrapolateBC(a_state, valid, a_dx, i, Side::Lo, a_order);
              break;
              
            case 3:
              DiriBC(a_state, valid, a_dx, a_homogeneous, PoissonBCLinearAlongZ, i, Side::Lo);
              break;
              
            default:
              MayDay::Error("bogus bc flag lo");
              break;
          }
        
        if (!a_domain.domainBox().contains(ghostBoxHi))
          
          switch (bcHi[i]) {
            case 0:
              DiriBC(a_state, valid, a_dx, a_homogeneous, PoissonBCParseValue, i, Side::Hi);
              break;
              
            case 1:
              NeumBC(a_state, valid, a_dx, a_homogeneous, PoissonBCParseValue, i, Side::Hi);
              break;
              
            case 2:
              ExtrapolateBC(a_state, valid, a_dx, i, Side::Hi, a_order);
              break;
              
            case 3:
              DiriBC(a_state, valid, a_dx, a_homogeneous, PoissonBCLinearAlongZ, i, Side::Hi);
              break;
              
            default:
              MayDay::Error("bogus bc flag hi");
              break;
          }
        
      }
    // end if is not periodic in ith direction
  }
}

void field::define(DisjointBoxLayout a_grids, IntVect ivGhost1, IntVect ivGhost2) {
  ParmParse pp;
  pp.get("verbosity",  s_verbosity);
  m_E.define(a_grids, SpaceDim, ivGhost1);
  m_Emag.define(a_grids, 1, ivGhost1);
  m_EEdge.define(a_grids, 1, ivGhost2);
}

void field::define(const field& a_field) {
  DisjointBoxLayout grids = a_field.m_Emag.disjointBoxLayout();
  IntVect ivGhost1, ivGhost2;
  ivGhost1 = a_field.m_Emag.ghostVect();
  ivGhost2 = a_field.m_EEdge.ghostVect();
  define(grids, ivGhost1, ivGhost2);
}

void field::copyTo(field& des) {
  m_EEdge.copyTo(m_EEdge.interval(), des.m_EEdge, des.m_EEdge.interval());
  m_E.copyTo(m_E.interval(), des.m_E, des.m_E.interval());
  m_Emag.copyTo(m_Emag.interval(), des.m_Emag, des.m_Emag.interval());
}

#include "NamespaceFooter.H"
