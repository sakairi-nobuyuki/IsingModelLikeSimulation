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


    //  allocate spin and spin record book
    s.resize (n_particles);
    s_record_book.resize (n_particles);
    

    //  initializing random spin
    initRandomSpins (); 

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
            
            Peval = 0.0;
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



void parallelArray1D::obtainCyclicBoundaryCondition () {
    //  cyclic boundary condition 
    s[0] = s[n_particles - 2];
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