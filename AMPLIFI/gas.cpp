//
//  gas.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#include <algorithm>
#include "gas.hpp"

using namespace std;

gas::gas(std::string a_name, bool a_uniformity, Real a_N, int a_numOfIonspe, double diffCoef)
: m_name(a_name), m_N(a_N), m_uniformity(a_uniformity), m_numOfIonSpe(a_numOfIonspe), m_elecDiffCoef(diffCoef)
{}

gas::gas(const gas& a_gas)
: m_name(a_gas.m_name), m_N(a_gas.m_N), m_uniformity(a_gas.m_uniformity), m_numOfIonSpe(a_gas.m_numOfIonSpe), m_elecDiffCoef(a_gas.m_elecDiffCoef) {

  EDpdentProcs = a_gas.EDpdentProcs;
  
  if(a_gas.bgdDensityProfile)
    //bgdDensityProfile = new parameterizedFunction(*a_gas.bgdDensityProfile);
    bgdDensityProfile = a_gas.bgdDensityProfile->clone();
  else
    bgdDensityProfile = NULL;
  
  if(a_gas.densityFileIF)
    densityFileIF = new DataFileIFReduced(*a_gas.densityFileIF);
  else
    densityFileIF = NULL;
}

gas& gas::operator=(const gas& a_gas) {
  m_name = a_gas.m_name;
  m_N = a_gas.m_N;
  m_uniformity = a_gas.m_uniformity;
  m_numOfIonSpe = a_gas.m_numOfIonSpe;
  m_elecDiffCoef = a_gas.m_elecDiffCoef;
  EDpdentProcs = a_gas.EDpdentProcs;
  
  if(a_gas.bgdDensityProfile)
    //bgdDensityProfile = new parameterizedFunction(*a_gas.bgdDensityProfile);
    bgdDensityProfile = a_gas.bgdDensityProfile->clone();
  else
    bgdDensityProfile = NULL;
  
  if(a_gas.densityFileIF)
    densityFileIF = new DataFileIFReduced(*a_gas.densityFileIF);
  else
    densityFileIF = NULL;
  
  return *this;
}

gas::~gas() {
  
  if (bgdDensityProfile != NULL)
    delete bgdDensityProfile;
  
  if (densityFileIF != NULL)
    delete densityFileIF;
}
  
double gas::getBackgroundDensity(vector<double> point) {
  
  if (bgdDensityProfile != NULL)
    return bgdDensityProfile->value(point);
  else if (densityFileIF != NULL)
    return densityFileIF->value(point);
  else {
    cerr << "gas background density missing!" << endl;
    abort();
  }
}
