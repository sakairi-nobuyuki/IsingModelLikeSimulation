#include <iostream>
#include <math.h>
#include "2d_parallel_array.h"


int main () {
    parallelArray2D Domain1;
    

    //cout << "initialize Hamiltonian" << endl;
    Domain1.obtainHamiltonian ();
    //cout << "finished to initialize Hamiltonian" << endl;

    
    Domain1.printHamiltonianAndSpinStatus ();

    //  obtain distribution function Z


    //  obtain probability p from H and Z



    //  perturvation against H with Gibbs sampling


    return 0;
}