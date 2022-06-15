//
//  gas.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef gas_hpp
#define gas_hpp

// This is a simple class for representing gases. A key member is the class to
// describe a process like ionization, which is a piece function of E/N.

#include <string>
#include <vector>
#include <map>

#include "dataFileIFReduced.hpp"
#include "parameterizedFunction.hpp"
#include "NamespaceHeader.H"

class gas;
class processCoefficients;


double getProcessCoeff(double E, gas& g, std::string procName);


// The following is a list of pre-defined functions

// A*exp(-B/E)
double piecewiseExponential(double E, processCoefficients& coeff);
// A + B*E
double piecewiseLinear(double E, processCoefficients& coeff);
// A + B/E
double piecewiseReciprocalLinear(double E, processCoefficients& coeff);

class processCoefficients {
  
private:
  
public:
  int numOfPieces;
  int numOfCoeffs; // number of coefficients for each piece
  vector<double> Xlb; // the left bound of the interval of each piece
  vector<vector<double>> A; // the coefficients
    
  processCoefficients() {};
  processCoefficients(vector<double> xlb, vector<vector<double>> a);
  processCoefficients& operator=(const processCoefficients&);
  void define(vector<double> xlb, vector<vector<double>> a);
  void define(const processCoefficients &);
  virtual ~processCoefficients() {};
};


typedef double (*processFunc) (double EN, processCoefficients&);

class process {
  
private:
  
public:
  processCoefficients coeff;          // the coefficients for the function
  processFunc func = NULL;            // the function
  
  process() {};
  process(const processCoefficients&, const processFunc&);
  process& operator=(const process&);
};


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
  
  std::map<std::string, process> processes;    // the name and process pair
  parameterizedFunction*         bgdDensityProfile = NULL;
  DataFileIFReduced*             densityFileIF = NULL;
  std::map<std::string, piecewiseFunction> EDpdentProcs;
};

#include "NamespaceFooter.H"

#endif /* gas_hpp */

