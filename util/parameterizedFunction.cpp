//
//  parameterizedFunction.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 6/10/22.
//

#include <cmath>
#include <iostream>
#include <cfloat>
#include <algorithm>

#include "parameterizedFunction.hpp"

map<string, functionPointer> singleFunction::functionCollection = {
  {"const",         constantValue},
  {"exp",           exponential},
  {"rcpExp",        reciprocalExponential},
  {"linear",        linear},
  {"rcpLinear",     reciprocalLinear},
  {"stdAtm",        standardAtmosphere}
};

double constantValue(vector<double>& x, vector<double>& p) {
  return p[0];
}

// e.g., N0*exp((z-z0)/H)
double exponential(vector<double>& x, vector<double>& p) {
  if (fabs(p[2]) > DBL_MIN)
    return p[0]*exp((x.back()-p[1])/p[2]);
  else {
    cerr << "division by zero in exponential()!" << endl;
    return 0.0;
  }
}

// e.g., A*exp(B/E)
double reciprocalExponential (vector<double>& x, vector<double>& p) {
  if (fabs(x.back()) > DBL_MIN)
    return p[0]*exp(p[1]/x.back());
  else {
    cerr << "division by zero in reciprocalExponential()!" << endl;
    return 0.0;
  }
}

// e.g., A+B*E
double linear (vector<double>& x, vector<double>& p) {
  return p[0] + p[1]*x.back();
}

// e.g., A+B/E
double reciprocalLinear (vector<double>& x, vector<double>& p) {
  if (fabs(x.back()) > DBL_MIN)
    return p[0] + p[1]/x.back();
  else {
    cerr << "division by zero reciprocalLinear!" << endl;
    return 0.0;
  }
}

double Satm(double h);
double standardAtmosphere(vector<double>& x, vector<double>& p) {
  double h = (x.back()*p[1] - p[0])/1000;
  return Satm(h)*2.688e25/p[2];
}

singleFunction::singleFunction(const string& name, const vector<double>& p)
: funcName(name), funcParam(p) {
  func = functionCollection[funcName];
}

singleFunction::singleFunction(const functionPointer& f, const vector<double>& p)
: func(f), funcParam(p) {}

singleFunction::singleFunction(const string& name, const functionPointer& f, const vector<double>& p)
: funcName(name), func(f), funcParam(p) {}

singleFunction::singleFunction(const singleFunction& f)
: funcName(f.funcName), func(f.func), funcParam(f.funcParam) {}

singleFunction& singleFunction::operator=(const singleFunction& f) {
  funcName  = f.funcName;
  funcParam = f.funcParam;
  func      = f.func;
  return *this;
}

double singleFunction::value(double x) {
  std::vector<double> xvec = {x};
  return value(xvec);
}

double singleFunction::value(vector<double> x) {
  return func(x, funcParam);
}


piecewiseFunction& piecewiseFunction::operator=(const piecewiseFunction& f) {
  numPieces  = f.numPieces;
  leftBound  = f.leftBound;
  funcs      = f.funcs;
  
  return *this;
}

piecewiseFunction::piecewiseFunction(const int num, const vector<double>& lb, const vector<string>& fNames, const vector<vector<double>>& p) : numPieces(num), leftBound(lb) {
  funcs.resize(num);
  for(int i = 0; i != numPieces; i++)
    funcs[i] = singleFunction(fNames[i], p[i]);
}

double piecewiseFunction::value(double x) {
  std::vector<double> xvec = {x};
  return value(xvec);
}

double piecewiseFunction::value(vector<double> x) {
  double eps = 1e-20;
  double xp = x.back();
  
  xp = xp + eps;
  vector<double>::iterator it = find_if(leftBound.begin(), leftBound.end(), greaterOrEqual(xp));
  it--;
  long index = it-leftBound.begin();
  
  if (index > -1)
    return funcs[index].value(x);
  else {
    cerr << "the input is out of the valid range of the piecewise function !" << endl;
    return 0.0;
    //abort();
  }
}

double Satm(double h) {
  double result, dz, nhi, nlo;
  int hi, lo, med;
  
  /* US Standard Atmosphere:
    We replaced original 2.5e19cm-3 at 0 km altitude
    by 2.688e19 cm-3 which is our reference
    number density at temperature 273 K.
    */
  
  double zz[]={0e+0, 5e+0,  1e+1, 1.5e+1,        2e+1,
       2.5e+1, 3e+1,  3.5e+1,   4e+1,   4.5e+1,
       5e+1,  5.5e+1,   6e+1,   6.4e+1, 6.8e+1,
       7.2e+1,   7.6e+1, 8e+1,  8.4e+1,   8.8e+1,
       9.2e+1, 9.6e+1, 1e+2,  1.08e+2, 1.14e+2,
       1.2e+2, 1.26e+2, 1.32e+2,  1.4e+2,   1.5e+2};
  
  
  double n[]={2.688e+19, 1.53e+19, 8.59e+18,4.05e+18,1.85e+18,
    8.33e+17, 3.83e+17,1.76e+17,8.31e+16,4.088e+16,
    2.13e+16, 1.181e+16, 6.33e+15, 3.93e+15, 2.39e+15,
    1.39e+15,7.72e+14,4.03e+14,1.99e+14, 9.48e+13,
    4.37e+13, 2.07e+13,1.04e+13,3.18e+12,1.43e+12,
       6.61e+11, 3.4e+11,1.91e+11,9.7e+10, 4.92e+10};
  
  if (h<0 || h>150) cerr << "h is outside of 0-150 km in standard atmosphere";
  
  /* Find interval to which given h belongs */
  lo=0;
  hi=29;
  while(hi>lo+1) {
    med=(hi+lo)/2;
    dz=h-zz[med];
    if(dz<0)  hi=med;
    if(dz>=0) lo=med;
  }
  
  nhi=log10(n[hi]);
  nlo=log10(n[lo]);
  
  result=nlo+(nhi-nlo)*(h-zz[lo])/(zz[hi]-zz[lo]);
  result=pow(10.,result)/n[0];
  return result;
}


