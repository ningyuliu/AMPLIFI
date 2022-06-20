//
//  parameterizedFunction.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 6/10/22.
//

#ifndef parameterizedFunction_hpp
#define parameterizedFunction_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <map>

using namespace std;

double constantValue(vector<double>& x, vector<double>& p);

// y = p[0] * exp( (x.back() - p[1]) / p[2] ); e.g., N0*exp((z-z0)/H)
double exponential (vector<double>& x, vector<double>& p);

// y = p[0] * exp( p[1] / x.back() ); e.g., A*exp(-B/E)
double reciprocalExponential (vector<double>& x, vector<double>& p);

// y = p[0] + p[1] * x.back(); e.g., A + B*E
double linear (vector<double>& x, vector<double>& p);

// y = p[0] + p[1] / x.back(); e.g., A + B/E
double reciprocalLinear (vector<double>& x, vector<double>& p);

// standard atmosphere
double standardAtmosphere (vector<double>& x, vector<double>& p);

// Wait ionosphere profile
double ionosphereWait(vector<double>& x, vector<double>& p);

double ionosphereTanh(vector<double>& x, vector<double>& p);

typedef double (*functionPointer) (vector<double>& x, vector<double>& parameter);

class parameterizedFunction {
  
private:
  
public:
  virtual ~parameterizedFunction() {};
  virtual parameterizedFunction* clone() const = 0;
  virtual double value(vector<double> x) = 0;
  
};

class singleFunction : public parameterizedFunction {
  
private:
  
public:
  singleFunction() {};
  singleFunction(const string& name, const vector<double>& p);
  singleFunction(const functionPointer& f, const vector<double>& p);
  singleFunction(const string& name, const functionPointer& fp, const vector<double>& p);
  singleFunction(const singleFunction& pf);
  singleFunction& operator=(const singleFunction&);
  virtual ~singleFunction() {};
  
  virtual singleFunction* clone() const {
    return new singleFunction(*this);
   }
  
  double value(vector<double> x);
  double value(double x);
  
  string              funcName;
  functionPointer     func = NULL;        // the function
  vector<double>      funcParam;          // the parameters for the function
  
  static map<string, functionPointer> functionCollection;
};



// A Functor
class greaterOrEqual {
  
private:
  double m_x;
  
public:
  greaterOrEqual(double a_x) : m_x(a_x) {  }
  
  // This operator overloading enables calling operator function () on objects of comparison
  int operator () (const double xlb) const {
    return xlb >= m_x;
  }
};



class piecewiseFunction : public parameterizedFunction {
  
private:
  
public:
  piecewiseFunction() {};
  piecewiseFunction(const int num, const vector<double>& lb, const vector<string>& funcNames,
                    const vector<vector<double>>& paramVect);
  piecewiseFunction& operator=(const piecewiseFunction&);
  virtual ~piecewiseFunction() {};
  
  virtual piecewiseFunction* clone() const {
    return new piecewiseFunction(*this);
   }
  
  int                             numPieces;
  vector<double>                  leftBound;
  vector<singleFunction>          funcs;

  double value(double x);
  double value(vector<double> x);
};




#endif /* parameterizedFunction_hpp */
