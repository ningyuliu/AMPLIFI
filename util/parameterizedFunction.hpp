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

// y = p[0] * exp( (p[1]-x.last()) / p[2] )
double exponential (vector<double>& x, vector<double>& p);

typedef double (*functionPointer) (vector<double>& x, vector<double>& parameter);


class parameterizedFunction {
  
private:
  
public:
  parameterizedFunction() {};
  parameterizedFunction(const string& name, const vector<double>& p);
  parameterizedFunction(const functionPointer& f, const vector<double>& p);
  parameterizedFunction(const string& name, const functionPointer& fp, const vector<double>& p);
  parameterizedFunction(const parameterizedFunction& pf);
  parameterizedFunction& operator=(const parameterizedFunction&);
  virtual ~parameterizedFunction() {};
  
  double value(vector<double> x);
  
  string              funcName;
  functionPointer     func = NULL;        // the function
  vector<double>      funcParam;          // the parameters for the function
  
  static map<string, functionPointer> functionCollection;
};

#endif /* parameterizedFunction_hpp */
