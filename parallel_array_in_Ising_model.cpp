#include <iostream>
#include <math.h>
#include "2d_parallel_array.h"
#include "1d_parallel_array.h"


int main () {
    int i;
    parallelArray1D Domain1;
    
    //  perturvation against H with Gibbs sampling

    Domain1.executeGibbsSampling ();


    return 0;
}