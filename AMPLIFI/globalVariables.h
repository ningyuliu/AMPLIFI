//
//  globalVariables.h
//  AMPLIFI
//
//  Created by Ningyu Liu on 4/28/22.
//

#ifndef globalVariables_h
#define globalVariables_h

// global constants for normalization

namespace normalization {

  // this is meant to minimize the number of changes required in the inputs file
  // for simulation at different altitudes; input quantities are scaled
  // according to similarity laws using this
  extern double scalingFactor;
  
  extern double EBar;
  extern double muBar;
  extern double lBar;
  extern double tBar;
  extern double nBar;
  extern double phiBar;

}

namespace numerical {
  
// In principle, Ac doesn't have to a global, as the cfl number is passed as an
// to generate the AMRLevelAdvecDiffuse class.
  extern double Ac;
  extern double Ai;
  extern double Ad;

}

#endif /* globalVariables_h */
