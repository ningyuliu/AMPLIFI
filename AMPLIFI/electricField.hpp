//
//  electricField.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef electricField_hpp
#define electricField_hpp

// This is for representing electric field and related boundary conditions

#include <functional>
#include <cmath>

#include "RealVect.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "BCFunc.H"

#include "NamespaceHeader.H"

// define a BC with constant value
extern void PoissonBCParseValueDiri(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);

extern void PoissonBCParseValueNeum(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);

// const potential at low and high of the last dimension, but linearly varies
// along the last dimension on other boundaries, where bc_lo/hi stores
// the electric field (- changing rate of the potential).
extern void PoissonBCLinearAlongZ(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values);

extern void EPotParseBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous);

extern Real linearTimeLinearZ(Real* pos, int* dir, Side::LoHiSide* side, Real time);

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


class TimeBCValueFunction : public BCValueFunction
{
public:
    using ValueFunc = std::function<Real(Real* pt, int* dir, Side::LoHiSide* side, Real time)>;

    TimeBCValueFunction() : m_time(0.0)
    {
        // Default behavior: sine wave in time
        m_func = [](Real* pt, int* dir, Side::LoHiSide* side, Real time) {
            return sin(2.0 * M_PI * time);
        };
    }

    void setTime(Real a_time) { m_time = a_time; }

    void setValueFunction(ValueFunc func) { m_func = func; }

    void operator()(Real* a_pos, int* a_dir, Side::LoHiSide* a_side, Real* a_value) override
    {
        a_value[0] = m_func(a_pos, a_dir, a_side, m_time);
    }

    Real m_time;

private:
    ValueFunc m_func;
};

class TimeDependentBCFunction : public BCFunction
{
public:
    TimeDependentBCFunction()
    {
        m_valueFunc = RefCountedPtr<TimeBCValueFunction>(new TimeBCValueFunction());
        // use the holder to interface Chombo's existing BC implementations like DiriBC
        m_holder = BCValueHolder(m_valueFunc);
    }

    void operator()(FArrayBox&           a_state,
                    const Box&           a_valid,
                    const ProblemDomain& a_domain,
                    Real                 a_dx,
                    bool                 a_homogeneous) override;

    void setTime(Real a_time)
    {
      m_valueFunc->setTime(a_time);
      pout() << "calling TimeDependentBCFunction.setTime() " << endl;
    }

    void setValueFunction(TimeBCValueFunction::ValueFunc func)
    {
        m_valueFunc->setValueFunction(func);
    }

private:
    RefCountedPtr<TimeBCValueFunction> m_valueFunc;
    BCValueHolder m_holder;
};

inline void setTimeHelper(BCHolder& bcHolder, Real time)
{
    RefCountedPtr<BCFunction> bcfunc = bcHolder.getBCFunction();
    if (bcfunc != nullptr)
    {
        // Use the conversion operator to get the raw pointer
      TimeDependentBCFunction* timeDepBC = dynamic_cast<TimeDependentBCFunction*>(bcfunc.operator->());

        if (timeDepBC != nullptr)
        {
          timeDepBC->setTime(time);
          pout() << "inside setTimeHelper" << endl;
        }
    }
}

#include "NamespaceFooter.H"

#endif /* electricField_hpp */
