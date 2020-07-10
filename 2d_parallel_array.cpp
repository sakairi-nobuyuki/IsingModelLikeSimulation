#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "2d_parallel_array.h"

using namespace std;

parallelArray2D::parallelArray2D () {
    n_spin_state = 3;
    n_x   = 100;
    n_y   = 2;
    n_paticles = n_x * n_y;

    //  allocate spin and spin record book
    s.resize (n_paticles);
    s_record_book.resize (n_paticles);
    Sigma.resize (n_paticles);

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
    for (i = 0; i < n_paticles; i++)  s[i] = rand () % n_spin_state;
}


void parallelArray2D::printSpinState () {
    cout << "state of spin" << endl;
    for (i_y = 0; i_y < n_y; i_y++) {
        cout << "  ";
        for (i_x = 0; i_x < n_x; i_x++) {   
            //cout << s[obtainOneDimPosFromTwoDimPos (i_x, i_y)];
            cout << Sigma[obtainOneDimPosFromTwoDimPos (i_x, i_y)];
        }
        cout << endl;

    }
}

void parallelArray2D::initSigma () {
    cout << "calc sigma" << endl;
    for (i_y = 0; i_y < n_y; i_y++) {
        for (i_x = 0; i_x < n_x; i_x++) {
            Sigma[obtainOneDimPosFromTwoDimPos (i_x, i_y)] 
                = obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x - 1, i_y)])
                + obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x + 1, i_y)])
                + obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x, i_y - 1)])
                + obtainSigmaFromSpin (s[obtainOneDimPosFromTwoDimPos (i_x, i_y + 1)]);
            cout << "  sigma " << i_x << " " << i_y  << " " << Sigma[obtainOneDimPosFromTwoDimPos (i_x - 1, i_y)] << endl;
        }
    }

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