#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <numeric>
#include "2d_parallel_array.h"

using namespace std;

parallelArray2D::parallelArray2D () {
    n_spin_state = 3;
    n_x   = 100;
    n_y   = 10;
    n_particles = n_x * n_y;

    //  allocate spin and spin record book
    s.resize (n_particles);
    s_record_book.resize (n_particles);
    Sigma.resize (n_particles);

    IonicStrength = 0.01;
    Phi   = 0.05;        //  volume fraction
    Delta = 1.0e-09;     //  width of the particle
    L     = 500.0e-09;   //  length of  th particle
    Cexcl = 0.8;

    kB = 1.0;
    T  = 100;
    Beta = 1.0 / (kB * T);
    
    //  initialize geometrical parameters
    Kappa = 0.304 / sqrt (IonicStrength) * 1.0e+09;
    //Dbar = Delta * (1.0 / Phi - 1.0);    ///  1 dimensional case
    Dbar = sqrt (2.0 * L * Delta / Phi) - 0.5 * Delta;
    Dtil = 0.5 * (L - Delta) * Delta / L * Cexcl;
    CvdW  = 0.0001;
    Alpha = 1.0e+09;
    
    VerPhiBar = exp (-Kappa * Dbar);
    VerPhiTil = exp (-Kappa * Dtil);
    VerPhi    = VerPhiBar * (pow (VerPhiTil, 2) - VerPhiTil);
    

    cout << "initializin geometrical parameters" << endl;
    cout << "  phi:   " << Phi << "(-)" << endl;
    cout << "  kappa: " << Kappa * 1.0e+09 << "(m+9)" << endl;
    cout << "  delta: " << Delta * 1.0e+09 << "(m-9)" << endl;
    cout << "  l:     " << L * 1.0e+09 << "(m-9)" << endl;
    cout << "  Dbar:  " << Dbar * 1.0e+09 << "(m-9)" << endl;
    cout << "  Dtil:  " << Dtil * 1.0e+09 << "(m-9)" << endl;
    cout << "  Kappa Dbar:    " << Kappa * Dbar << "(-)" << endl;
    cout << "  Kappa Dtil:    " << Kappa * Dtil << "(-)" << endl;
    cout << "  VerPhiBar:     " << VerPhiBar << "(-)" << endl;
    cout << "  VerPhiTil:     " << VerPhiTil << "(-)" << endl;
    cout << "  IonicStrength: " << IonicStrength << "(mol/L)" << endl;
    cout << "  Alpha:   " << Alpha << "(-)" << endl << endl;

    // initialization with random spins
    initRandomSpins ();
    initSigma ();


}


void parallelArray2D::initRandomSpins () {
    for (i = 0; i < n_particles; i++)  s[i] = rand () % n_spin_state;
}


void parallelArray2D::printSpinState () {
    cout << "state of spin" << endl;
    for (i_y = 0; i_y < n_y; i_y++) {
        cout << "  ";
        for (i_x = 0; i_x < n_x; i_x++) {   
            cout << s[obtainOneDimPosFromTwoDimPos (i_x, i_y)];
            //cout << Sigma[obtainOneDimPosFromTwoDimPos (i_x, i_y)];
        }
        cout << endl;

    }
}

void parallelArray2D::initSigma () {
    //cout << "calc sigma" << endl;
    for (i_y = 0; i_y < n_y; i_y++) {
        for (i_x = 0; i_x < n_x; i_x++) {
            Sigma[obtainOneDimPosFromTwoDimPos (i_x, i_y)] 
                = obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x - 1, i_y)])
                + obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x + 1, i_y)])
                + obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x, i_y - 1)])
                + obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x, i_y + 1)]);
        //    cout << "  sigma " << i_x << " " << i_y  << " " << Sigma[obtainOneDimPosFromTwoDimPos (i_x - 1, i_y)] << endl;
        }
    }

}

void  parallelArray2D::executeGibbsSampling2D () {
    int i_sample, n_shifted_particles;
    double Prand, Peval;

    //  initialize spin recordbook
    cout << "start perturbation" << endl << "initializing spin book" << endl;
    for (i = 0; i < n_particles; i++) s_record_book[i] = 0;

    //  choose particles
    for (i = 0; n_shifted_particles = accumulate (s_record_book.begin (), s_record_book.end (), 0) < n_particles; i++) {
        i_sample = rand () % n_particles;
        //cout << i_sample << "shall be evaluated first" << endl;
        if (s_record_book[i_sample] == 0) {
            
            Peval = 0.0;
            Prand = (double) (rand () % 100) * 0.01;
            cout << i << "th sampling: " << endl;
            cout << "  " << i_sample << " th spin, c = " << s[i_sample] << endl;
            cout << "  Pcond = " << Peval << ", Prand = " << Prand << endl;
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

void parallelArray2D::obtainJointProbability () {
    obtainHamiltonian ();
    cout << "joint probability: " << exp (-Beta * H);

}

void parallelArray2D::obtainHamiltonian () {
    //for (i = 0; i < n_paticles; i++) {
    cout << "going to obtain Hamiltonian" << endl;
    applyCyclicBoundaryCondition ();
    cout << "finished cyclic Boundary Condition" << endl;
    H = 0.0;
    for (i = 0; i < n_particles - (n_x + 1); i++) {        
        //cout << i << endl;
        H += (exp (-Kappa * Dtil * (obtainSigmaFromSpin (s[i]) + obtainSigmaFromSpin (s[i + 1])))
            + exp (-Kappa * Dtil * (obtainSigmaFromSpin (s[i]) + obtainSigmaFromSpin (s[i + n_x]))))
            * (i % n_x != 0);
    }
    H *= exp (-Kappa * Dbar);
    cout << "Hamiltonian: " << H << endl;
}

void parallelArray2D::applyCyclicBoundaryCondition () {
    cout << "cyclic boundary condition" << endl;
    for (i_x = 0; i_x < n_x - 1; i_x++)  {
        //cout << i_x << " " << obtainOneDimPosFromTwoDimPos (i_x, n_x - 2) << " " << n_x << " " << n_y << endl;
        s[obtainOneDimPosFromTwoDimPos (i_x, 0)]       = s[obtainOneDimPosFromTwoDimPos (i_x, n_y - 2)]; 
    }
    
    for (i_x = 0; i_x < n_x - 1; i_x++)  {
        s[obtainOneDimPosFromTwoDimPos (i_x, n_y - 1)] = s[obtainOneDimPosFromTwoDimPos (i_x, 1)]; 
    }
    for (i_y = 0; i_y < n_y - 1; i_y++)  {
        s[obtainOneDimPosFromTwoDimPos (0, i_y)]       = s[obtainOneDimPosFromTwoDimPos (n_x - 2, i_y)]; 
    }
    for (i_y = 0; i_y < n_y - 1; i_y++)  {
        s[obtainOneDimPosFromTwoDimPos (n_x - 1, i_y)] = s[obtainOneDimPosFromTwoDimPos (1, i_y)];         
    }
    printSpinState ();
}


inline double parallelArray2D::obtainCondProb23 (int i) {
    return 1.0 / (2.0 + exp (Beta * obtainSigmaFromSpin (s[i]) * VerPhiBar * (VerPhi - pow (VerPhiTil, 2))));
}

inline double parallelArray2D::obtainCondProb1 (int i) {
    return 1.0 / (1.0 + 2.0 * exp(Beta * obtainSigmaFromSpin (s[i]) * VerPhiBar * (pow (VerPhiTil, 2) - VerPhiTil)));
}

inline double parallelArray2D::obtainSigmaFromSpin (double s) {
    return VerPhiBar * VerPhiTil * ((s != 2) * 1.0 + (s == 2) * VerPhiTil);
}

inline int parallelArray2D::obtainOneDimPosFromTwoDimPos (int i_x, int i_y) {
    return i_y * n_x + i_x;

}