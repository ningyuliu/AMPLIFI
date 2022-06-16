//
//  gas.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef gas_hpp
#define gas_hpp

// This is a simple class for representing gases.

#include <string>
#include <vector>
#include <map>

#include "dataFileIFReduced.hpp"
#include "parameterizedFunction.hpp"

class gas {
  
private:
  
public:
  
  gas() {};
  gas(std::string a_name, bool uniformity, Real a_N, int a_numOfIonspe=1, double diffCoef=0);
  gas(const gas&);
  gas& operator=(const gas&);
  virtual ~gas();
  
  double getBackgroundDensity(vector<double> point);
  
  std::string             m_name;
  double                  m_N;                 // average/reference density
  bool                    m_uniformity;        // uniform density?
  int                     m_numOfIonSpe;       // electrons excluded
  double                  m_elecDiffCoef;
  
  parameterizedFunction*         bgdDensityProfile = NULL;
  DataFileIFReduced*             densityFileIF = NULL;
  std::map<std::string, piecewiseFunction> EDpdentProcs;
};

#endif /* gas_hpp */

