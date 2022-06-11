//
//  parameterizedFunction.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 6/10/22.
//

#include <math.h>
#include "parameterizedFunction.hpp"

map<string, functionPointer> parameterizedFunction::functionCollection = {
  {"exp", exponential}
};

parameterizedFunction::parameterizedFunction(const string& name, const vector<double>& p)
: funcName(name), funcParam(p) {
  func = functionCollection[funcName];
}

parameterizedFunction::parameterizedFunction(const functionPointer& f, const vector<double>& p)
: func(f), funcParam(p) {}

parameterizedFunction::parameterizedFunction(const string& name, const functionPointer& f, const vector<double>& p)
: funcName(name), func(f), funcParam(p) {}

parameterizedFunction::parameterizedFunction(const parameterizedFunction& pf)
: funcName(pf.funcName), func(pf.func), funcParam(pf.funcParam) {}

parameterizedFunction& parameterizedFunction::operator=(const parameterizedFunction& pf) {
  funcName  = pf.funcName;
  funcParam = pf.funcParam;
  func      = pf.func;  
  return *this;
}

double parameterizedFunction::value(vector<double> x) {
  return func(x, funcParam);
}

double exponential(vector<double>& x, vector<double>& p) {
  double y;
  y = p[0]*exp((x.back()-p[1])/p[2]);
  return y;
}
