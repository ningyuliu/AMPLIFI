//
//  gas.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#include <algorithm>
#include "gas.hpp"

#include "NamespaceHeader.H"

using namespace std;

double getProcessCoeff(double E, gas& g, string name) {
  return g.processes.at(name).func(E, g.processes.at(name).coeff);
}

processCoefficients::processCoefficients(vector<double> xlb, vector<vector<double>> a) : Xlb(xlb), A(a) {
  numOfPieces = Xlb.size();
  numOfCoeffs = A.size();
}

void processCoefficients::define(vector<double> xlb, vector<vector<double>> a) {
  Xlb = xlb;
  A = a;
  numOfPieces = Xlb.size();
  numOfCoeffs = A.size();
}

void processCoefficients::define(const processCoefficients & procCoeff) {
  define(procCoeff.Xlb, procCoeff.A);
}

processCoefficients& processCoefficients::operator=(const processCoefficients& procCoeff) {
  define(procCoeff);
  return *this;
}

process::process(const processCoefficients& c, const processFunc& f) {
  coeff = c;
  func = f;
}

process& process::operator=(const process& p) {
  coeff = p.coeff;
  func = p.func;
  
  return *this;
}


gas::gas(std::string a_name, bool a_uniformity, Real a_N, int a_numOfIonspe)
: m_name(a_name), m_N(a_N), m_uniformity(a_uniformity), m_numOfIonSpe(a_numOfIonspe)
{}

gas::gas(const gas& a_gas)
: m_name(a_gas.m_name), m_N(a_gas.m_N), m_uniformity(a_gas.m_uniformity), m_numOfIonSpe(a_gas.m_numOfIonSpe){
  
  processes = a_gas.processes;
  
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
  processes = a_gas.processes;
  
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

void gas::define (std::string a_name, bool a_uniformity, Real a_N, int a_numOfIonSpe) {
  m_name = a_name;
  m_uniformity = a_uniformity;
  m_N = a_N;
  m_numOfIonSpe = a_numOfIonSpe;
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

//#if __cplusplus > 199711L
//
//// implementation of A*exp(-B/E) using lambda expresssion
//double piecewiseExponential(const double& E, const processCoefficients& coeff) {
//  double eps = 1e-10;
//  double Emod = E + eps;
//  // use lambda expression to capture E
//  vector<const double>::iterator it = find_if(coeff.Xlb.begin(), coeff.Xlb.end(), [E](double xlb){ return xlb >= E; });
//  it--;
//
//  int index = it-coeff.Xlb.begin();
//  double y;
//  y = *(coeff.A[0].begin()+index) * exp(-*(coeff.A[1].begin()+index)/Emod);
//
//  return y;
//}
//
//// implementation of A + B/E
//double piecewiseLinear(const double& E, const processCoefficients& coeff) {
//  double eps = 1e-10;
//  double Emod = E + eps;
//  // use lambda expression to capture E
//  vector<const double>::iterator it = find_if(coeff.Xlb.begin(), coeff.Xlb.end(), [E](double xlb){ return xlb >= E; });
//  it--;
//
//  int index = it-coeff.Xlb.begin();
//  double y;
//  y = *(coeff.A[0].begin()+index) + *(coeff.A[1].begin()+index)/Emod;
//
//  return y;
//}
//
//#else

// A Functor
class compareToNumber
{
private:
  double m_x;
public:
  compareToNumber(double a_x) : m_x(a_x) {  }
  
  // This operator overloading enables calling operator function () on objects of comparison
  int operator () (const double xlb) const {
      return xlb >= m_x;
  }
};

// implementation of A*exp(-B/E)
double piecewiseExponential(double E, processCoefficients& coeff) {
  double eps = 1e-10;
  double Emod = E + eps;
  vector<double>::iterator it = find_if(coeff.Xlb.begin(), coeff.Xlb.end(), compareToNumber(E));
  it--;
  
  long index = it-coeff.Xlb.begin();
  double y;
  y = *(coeff.A[0].begin()+index) * exp(-*(coeff.A[1].begin()+index)/Emod);
  
  return y;
}

// implementation of A + B*E
double piecewiseLinear(double E, processCoefficients& coeff) {
  double eps = 1e-10;
  double Emod = E + eps;
  vector<double>::iterator it = find_if(coeff.Xlb.begin(), coeff.Xlb.end(), compareToNumber(E));
  it--;
  
  long index = it-coeff.Xlb.begin();
  double y;
  y = *(coeff.A[0].begin()+index) + *(coeff.A[1].begin()+index)*Emod;
  
  return y;
}

// implementation of A + B/E
double piecewiseReciprocalLinear(double E, processCoefficients& coeff) {
  double eps = 1e-10;
  double Emod = E + eps;
  vector<double>::iterator it = find_if(coeff.Xlb.begin(), coeff.Xlb.end(), compareToNumber(E));
  it--;
  
  long index = it-coeff.Xlb.begin();
  double y;
  y = *(coeff.A[0].begin()+index) + *(coeff.A[1].begin()+index)/Emod;
  
  return y;
}

//#endif

#include "NamespaceFooter.H"
