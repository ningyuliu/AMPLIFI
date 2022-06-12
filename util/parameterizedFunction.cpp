//
//  parameterizedFunction.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 6/10/22.
//

#include <math.h>
#include <iostream>

#include "parameterizedFunction.hpp"

map<string, functionPointer> singleFunction::functionCollection = {
  {"exp", exponential}
};

double exponential(vector<double>& x, vector<double>& p) {
  double y;
  y = p[0]*exp((x.back()-p[1])/p[2]);
  return y;
}

singleFunction::singleFunction(const string& name, const vector<double>& p)
: funcName(name), funcParam(p) {
  func = functionCollection[funcName];
}

singleFunction::singleFunction(const functionPointer& f, const vector<double>& p)
: func(f), funcParam(p) {}

singleFunction::singleFunction(const string& name, const functionPointer& f, const vector<double>& p)
: funcName(name), func(f), funcParam(p) {}

singleFunction::singleFunction(const singleFunction& pf)
: funcName(pf.funcName), func(pf.func), funcParam(pf.funcParam) {}

singleFunction& singleFunction::operator=(const singleFunction& pf) {
  funcName  = pf.funcName;
  funcParam = pf.funcParam;
  func      = pf.func;  
  return *this;
}

double singleFunction::value(vector<double> x) {
  return func(x, funcParam);
}



piecewiseFunction::piecewiseFunction(const int num, const vector<double>& lb, const vector<string>& fNames, const vector<vector<double>>& p) : numPieces(num), leftBound(lb) {
  
  paramFuncVect.resize(num);
  for(int i = 0; i != numPieces; i++)
    paramFuncVect[i] = singleFunction(fNames[i], p[i]);
}

double piecewiseFunction::value(vector<double> x) {
  
  double eps = 1e-20;
  double xp = x.back();
  
  xp = xp + eps;
  vector<double>::iterator it = find_if(leftBound.begin(), leftBound.end(), greaterOrEqual(xp));
  it--;
  long index = it-leftBound.begin();
  
  if (index > -1)
    return paramFuncVect[index].value(x);
  else {
    cerr << "the input is out of the valid range of the piecewise function !" << endl;
    return 0.0;
    //abort();
  }
}
