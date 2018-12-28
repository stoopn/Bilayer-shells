/*
 *  shell_io.h
 *  newshell_prj
 *
 *  Created by Norbert Stoop on 03.06.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SHELL_UTIL_H
#define SHELL_UTIL_H

#include "getpot.h"
#include "shellsystem.h"
#include "equation_systems.h"
#include "libmesh.h"
#include "mesh.h"
#include "o_string_stream.h"
#include <fstream>



//        "Polar" version without trigonometric calls

inline double randn_notrig(double mu=0.0, double sigma=1.0) {
  static bool deviateAvailable=false;                //        flag
  static float storedDeviate;                        //        deviate from previous calculation
  double polar, rsquared, var1, var2;

  //        If no deviate has been stored, the polar Box-Muller transformation is
  //        performed, producing two independent normally-distributed random
  //        deviates.  One is stored for the next round, and one is returned.

  if (!deviateAvailable) {
    //        choose pairs of uniformly distributed deviates, discarding those
    //        that don't fall within the unit circle

    do {
      var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
      var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
      rsquared=var1*var1+var2*var2;
    } while ( rsquared>=1.0 || rsquared == 0.0);
    //        calculate polar tranformation for each deviate
    polar=sqrt(-2.0*log(rsquared)/rsquared);
    //        store first deviate and set flag
    storedDeviate=var1*polar;
    deviateAvailable=true;
    //        return second deviate
    return var2*polar*sigma + mu;
  }
  //        If a deviate is available from a previous call to this function, it is
  //        returned, and the flag is set to false.
  else {
    deviateAvailable=false;
    return storedDeviate*sigma + mu;
  }
}

#endif
