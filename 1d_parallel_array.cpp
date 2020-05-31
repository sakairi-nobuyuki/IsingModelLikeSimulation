#ifndef __1D_PARALLEL_ARRAY__
#define __1D_PARALLEL_ARRAY__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <bitset>
#include <list>
#include <string>

#include "1d_parallel_array.h"
#include "utility.h"


using namespace std;

parallelArray1D::parallelArray1D () {
    n_status = 2;
    n_particles = 5;
    Kappa = 1.0;
    CvdW  = 0.0001;

    kB = 1.0;
    T  = 100;
    Beta = 1.0 / (kB * T);
    Z  = obtainDistributionFunction ();

    //  allocate spin and spin record book
    s.resize (n_particles);
    s_record_book.resize (n_particles);

    // initialization with random spins
    initRandomSpins ();

    Z = obtainDistributionFunction ();

}


double parallelArray1D::obtainDistributionFunction () {
    int i, j, k, spin_status_in_bin;
    string spin_state;
    list<string>::iterator all_spin_state_itr;
    
    setAllSpinState (all_spin_state);
    

    //  Z = sum exp (-Beta H (s))
    //cout << endl << endl << "all spinstate" << endl;
    for (all_spin_state_itr = all_spin_state.begin (); all_spin_state_itr != all_spin_state.end (); all_spin_state_itr++) {
        spin_state = *all_spin_state_itr;
        cout << spin_state << endl;
        
        for (i = 0; i < spin_state.size (); i++) {
            cout << "  " << spin_state[i] << endl;
        }
        
    }
        
    

    
    return 1.0;

}


void parallelArray1D::setAllSpinState (list<string> &AllSpinState) {
    int i;

    for (i = 0; i < pow (n_status, n_particles); i++) {
        string spin_status_in_str;
        stringstream bin_in_string;    

        bin_in_string << bitset <sizeof (i) * __CHAR_BIT__> (i);
        spin_status_in_str = bin_in_string.str ().substr (sizeof (i) * __CHAR_BIT__ - n_particles, n_particles);
        cout << i << " " << spin_status_in_str << endl;
        AllSpinState.push_back (spin_status_in_str);
    }

}


void parallelArray1D::initRandomSpins () {
    for (i = 0; i <n_particles; i++) {
        s[i] = rand () % n_status;
        cout << s[i] << " ";
    }
    cout << endl;
}

void parallelArray1D::perturbSpinDistributionWithGibbsSampling () {
    int i_sample, i_one_hot, j_one_hot;

    //  initialize spin recordbook
    for (i = 0; i <n_particles; i++) s_record_book[i] = 0;
        
    //  choose which spin to move
    

    //  evaluate if spin direction to be moved
    



}

void parallelArray1D::obtainHamiltonian (int *s, int n_particles) {
    // boundary condition
    cout << "set cyclic boundary condition" << endl;
    obtainCyclicBoundaryCondition ();
    printHamiltonianAndSpinStatus ();
    cout << "finished to set cyclic boundary condition" << endl;

    //  calculating Hamiltonian
    cout << "start to calculate Hamiltonian" << endl;
    H = 0.0;
    cout << "initialize Hamiltonian" << endl;
    for (i = 1; i < n_particles - 1; i++) {
        for (k = 0; k < n_particles; k++) {
            H += (obtainVanDerWaalsInteractionDiscrete1D (i, k) +
                obtainEDLInteractionDiscrete1D (i, k)) * s[i] * s[k];
        }
    }
}

double parallelArray1D::obtainProbability (double H) {
    return exp (-1.0 * Beta * H) / Z;
}

void parallelArray1D::obtainCyclicBoundaryCondition () {
    //  cyclic boundary condition 
    s[0] = s[n_particles - 2];
}


void parallelArray1D::printHamiltonianAndSpinStatus () {
    cout << "Hamiltonian: " << H << endl;
    cout << "Spin status:" << endl;
    for (i = 0; i <n_particles; i++) cout << s[i] << " ";
    cout << endl;
    cout << "End spin status:" << endl;
}



double parallelArray1D::obtainD1DistanceDiscrete1D (int x1, int x2) {
    if (x1 == x2) {
        return 0.0;
    } else {
        return 1.0;
    }
}


double parallelArray1D::obtainEDLInteractionDiscrete1D (int x1, int x2) {
    return obtainD1DistanceDiscrete1D (x1, x2) * exp (-1.0 * Kappa * obtainD1DistanceDiscrete1D (x1, x2));

}


double parallelArray1D::obtainVanDerWaalsInteractionDiscrete1D (int x1, int x2) {
    return -1.0 * obtainD1DistanceDiscrete1D (x1, x2) * CvdW / (obtainD1DistanceDiscrete1D (x1, x2) + 0.001);

}

double parallelArray1D::obtainD2DistanceDiscrete1D (int x1, int x2) {
    return fabs ((double) (x1 - x2));
}


#endif