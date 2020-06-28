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
#include <typeinfo>
#include <numeric>

#include "1d_parallel_array.h"
#include "utility.h"


using namespace std;

parallelArray1D::parallelArray1D () {
    n_status = 2;
    n_particles = 8;
    Kappa = 1.0;
    CvdW  = 0.0001;

    kB = 1.0;
    T  = 100;
    Beta = 1.0 / (kB * T);
    Z  = obtainDistributionFunction ();

    //  allocate spin and spin record book
    s.resize (n_particles);
    s_record_book.resize (n_particles);
    

    //  initializing random spin
    initRandomSpins (); 

}


double parallelArray1D::obtainDistributionFunction () {
    int i, j, k, spin_status_in_bin, spin_target[n_particles], spin_test;
    const char* spin_state_char;
    string spin_state, spin_state_1_digit;
    list<string> all_spin_state;
    list<string>::iterator all_spin_state_itr;
    vector<int> spin_subject (n_particles);
    double Z, H;

    Z = 0.0;
    H = 0.0;
    
    setAllSpinState (all_spin_state);
    
    //  Z = sum exp (-Beta H (s))
    //cout << endl << endl << "all spinstate" << endl;
    for (all_spin_state_itr = all_spin_state.begin (); all_spin_state_itr != all_spin_state.end (); all_spin_state_itr++) {
        spin_state = *all_spin_state_itr;    //  transform list of string to string
        //cout << "spin state: " << spin_state << ", size of spin state = " << spin_state.size () << endl;    
        for (i = 0; i < spin_state.size (); i++) {
            spin_state_1_digit = spin_state[i];
            //cout << spin_state_1_digit << " " << spin_state[i] << endl;
            spin_test = stoi (spin_state_1_digit);
            if (spin_test == 0) spin_test = -1;
            //cout << spin_test << endl;
            spin_subject[i] = spin_test;
            cout << spin_subject[i];
        }
        cout << endl;
        
        H = obtainHamiltonian (spin_subject);
        Z += exp (-Beta * H);

        cout << "H = " << H << endl;
        
    }
    
    cout << "Z = " << Z << endl;
    
    return Z;
}


void parallelArray1D::setAllSpinState (list<string>& AllSpinState) {
    int i;

    for (i = 0; i < pow (n_status, n_particles); i++) {
        string spin_status_in_str;
        stringstream bin_in_string;    

        bin_in_string << bitset <sizeof (i) * __CHAR_BIT__> (i);
        spin_status_in_str = bin_in_string.str ().substr (sizeof (i) * __CHAR_BIT__ - n_particles, n_particles);
        cout << i << " " << spin_status_in_str << endl;
        AllSpinState.push_back (spin_status_in_str);
    }
    cout << "all spin set finished" << endl;
}


void parallelArray1D::initRandomSpins () {
    for (i = 0; i <n_particles; i++) {
        s[i] = rand () % n_status;
        if (s[i] == 0)  s[i] = -1;
        //cout << s[i] << " ";
    }
    //cout << endl;
}

void parallelArray1D::perturbSpinDistributionWithGibbsSampling () {
    int i_sample, i_one_hot, j_one_hot, n_shifted_particles;
    double Peval, Prand;

    //  initialize spin recordbook
    cout << "start perturbation" << endl << "initializing spin book" << endl;
    for (i = 0; i <n_particles; i++) s_record_book[i] = 0;
        
    //  choose which spin to move
    for (i = 0; n_shifted_particles = accumulate (s_record_book.begin (), s_record_book.end (), 0) < n_particles; i++) {
        i_sample = rand () % n_particles;
        cout << i_sample << "shall be evaluated first" << endl;
        if (s_record_book[i_sample] == 0) {
            Peval = exp (-Beta * obtainHamiltonian (s)) / Z;
            Prand = (double) (rand () % 100) * 0.01;
            if (Peval < Prand) {
                s[i_sample] = 1;
                cout << "up " << Peval << " " << Prand << endl;
            } else {
                s[i_sample] = -1;
                cout << "down " << Peval << " " << Prand  << endl;
            }
            s_record_book[i_sample] = 1;
        } 

    }
    for (i = 0; i < s.size (); i++)  cout << s[i];
    cout << endl;
    //  evaluate if spin direction to be moved

}

double parallelArray1D::obtainHamiltonian (vector<int>& s) {
    int i;
    double H;

    // boundary condition
    //cout << "set cyclic boundary condition with spin state ";
    //for (i = 0; i < n_particles; i++) cout << s[i];
    //cout << endl;
    //cout << "finished to set cyclic boundary condition" << endl;

    //  calculating Hamiltonian
    //cout << "start to calculate Hamiltonian" << endl;
    H = 0.0;
    //cout << "initialize Hamiltonian" << endl;
    for (i = 0 ; i < s.size (); i++) {
        for (k = 0; k < i; k++) {
            J = obtainEDLInteractionDiscrete1D (i, k) * s[i] * s[k]
                - obtainVanDerWaalsInteractionDiscrete1D (i, k);
            H += J;
            cout << "i = " << i << ", j = " << k << ", s_i, s_k = " << s[i] << ", " << s[k] << ", J = " << J << ", p = " << exp (-Beta * J) << ", Z = " << Z << endl;
        }
    }
    
    
    return H;
}

double parallelArray1D::obtainProbability (double H) {
    return exp (-1.0 * Beta * H) / Z;
}

void parallelArray1D::obtainCyclicBoundaryCondition () {
    //  cyclic boundary condition 
    s[0] = s[n_particles - 2];
}


void parallelArray1D::printHamiltonianAndSpinStatus () {
    //cout << "Hamiltonian: " << H << endl;
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
    return obtainD1DistanceDiscrete1D (x1, x2) * exp (-1.0 * Kappa * obtainD2DistanceDiscrete1D (x1, x2));

}


double parallelArray1D::obtainVanDerWaalsInteractionDiscrete1D (int x1, int x2) {
    return -1.0 * obtainD1DistanceDiscrete1D (x1, x2) * CvdW / (obtainD2DistanceDiscrete1D (x1, x2) + 0.001);

}

double parallelArray1D::obtainD2DistanceDiscrete1D (int x1, int x2) {
    return fabs ((double) (x1 - x2));
}


#endif