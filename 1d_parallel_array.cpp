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
    n_status     = 2;
    n_particles  = 100;
    n_state_of_c = 10;
    
    IonicStrength = 0.01;
    Phi   = 0.1; 
    Delta = 1.0e-09;
    L     = 500.0e-09;
    CvdW  = 0.0001;
    Cex = 4.0;

    c_th  = 2;

    kB = 1.0;
    T  = 100;
    Beta = 1.0 / (kB * T);

    //  initialize geometrical parameters
    Kappa = 0.304 / sqrt (IonicStrength) * 1.0e-09;
    D   = Delta * (1.0 / Phi - 1.0);
    Dex = 0.5 * (L - Delta) * Delta / L * Cex;
    Alpha = 1.0;

    cout << "initializin geometrical parameters" << endl;
    cout << "  phi:   " << Phi << "(-)" << endl;
    cout << "  kappa: " << Kappa * 1.0e+09 << "(m-9)" << endl;
    cout << "  delta: " << Delta * 1.0e+09 << "(m-9)" << endl;
    cout << "  l:     " << L * 1.0e+09 << "(m-9)" << endl;
    cout << "  D:     " << D * 1.0e+09 << "(m-9)" << endl;
    cout << "  Dex:   " << Dex * 1.0e+09 << "(m-9)" << endl;
    cout << "  IonicStrength: " << IonicStrength << "(mol/L)" << endl;
    cout << "  Alpha:   " << Alpha << "(-)" << endl << endl;

    //  allocate spin and spin record book
    s.resize (n_particles);
    s_record_book.resize (n_particles);
    c.resize (n_particles);
    c_record_book.resize (n_particles);

    //  initializing random spin
    initRandomOrderParameter (); 

    //  transformation order parameter to spin
    transformOrderParameterToSpinState ();

    // confirmation
    cout << "Finished initialize order parameter and spin" << endl;
    printSpinAndOrderParameterState ();
    cout << endl;

}

void parallelArray1D::initRandomOrderParameter () {
    for (i = 0; i <n_particles; i++) {
        c[i] = rand () % n_state_of_c;
        //cout << c[i] << " ";
    }
    //cout << endl;
}


void parallelArray1D::printSpinAndOrderParameterState () {
    cout << "Order Parameter" << endl << "  ";
    for (i = 0; i < n_particles; i++) cout << c[i];
    cout << endl;

    cout << "Spin" << endl << "  ";
    for (i = 0; i < n_particles; i++) cout << s[i];
    cout << endl;

}

void parallelArray1D::executeGibbsSampling () {
    int i_sample, n_shifted_particles;
    double Prand, Peval;

    //  initialize spin recordbook
    cout << "start perturbation" << endl << "initializing spin book" << endl;
    for (i = 0; i <n_particles; i++) c_record_book[i] = 0;

    //  choose particles
    for (i = 0; n_shifted_particles = accumulate (c_record_book.begin (), c_record_book.end (), 0) < n_particles; i++) {
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
}

inline double parallelArray1D::obtainConditionalProb1D (int i_particle, int c_subject) {

    ProtoProb = 0.0;

    for (j = 0; j < n_state_of_c; j++) {
        J = obtainInteractionPotential (s[i_particle - 1], s[i_particle + 1]);
        ProtoProb += exp (-Beta * J * ((obtainSpinFromOrderParameter (c[i_particle]) - obtainSpinFromOrderParameter (c_subject)) 
            * (s[i_particle - 1] - s[i_particle + 1])));

    }

    return 1.0 / (1.0 );
}

inline double parallelArray1D::obtainInteractionPotential (int s_1, int s_2) {
    return -1.0 / D + Alpha * exp (-Kappa * (D - (s_1 + s_2) * Dex));
}

void parallelArray1D::transformOrderParameterToSpinState () {
    for (i = 0; i < n_particles; i++)   s[i] = ((c[i] - c_th) > 0);
}

inline int parallelArray1D::obtainSpinFromOrderParameter (int c) {
    return (c - c_th) > 0;
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