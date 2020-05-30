#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "2d_parallel_array.h"

using namespace std;

parallelArray2D::parallelArray2D () {

    LL = 20;
    Kappa = 1.0;
    CvdW  = 1.0;

    kB = 1.0;
    T  = 100;
    Beta = kB * T;
    Z  = obtainDistributionFunction ();

    //  allocate spin and spin record book
    s.resize (LL);
    for (i = 0; i < LL; i++)  s[i].resize (LL);
    s_record_book.resize (LL);
    for (i = 0; i < LL; i++) {
        s_record_book[i].resize (LL);
    }

    // initialization with random spins
    initRandomSpins ();



}

double parallelArray2D::obtainDistributionFunction () {
    //  Z = sum exp (-Beta H (s))
    return 1.0;
//    for (i = 0; i < LL; i++) {
//        for (j = 0; j < LL; j++) {
//            for (k = 0; k < LL; k++) {
//                for (l = 0; l < LL; l++) {

//                }
//            }
//        }

//    }
}


void parallelArray2D::initRandomSpins () {
    for (i = 0; i < LL; i++) {
        for (j = 0; j < LL; j++) {
            s[i][j] = rand () % 2;
            cout << s[i][j] << " ";
        }
        cout << endl;
    }

}

void parallelArray2D::perturbSpinDistributionWithGibbsSampling () {
    int i_sample, i_one_hot, j_one_hot;

    //  initialize spin recordbook
    for (i = 0; i < LL; i++) {
        for (j = 0; j < LL; j++)  s_record_book[i][j] = 0;
    }
    
    //  choose which spin to move
    for (i = 0; ; i++) {
        i_sample = rand () % (int) pow (LL, 2);
        i_one_hot = i_sample / LL;
        j_one_hot = i_sample % LL;
        if (s_record_book[i_one_hot][j_one_hot] == 1) continue;
        else break;
    }

    //  evaluate if spin direction to be moved




}

void parallelArray2D::obtainHamiltonian () {
    // boundary condition
    cout << "set cyclic boundary condition" << endl;
    obtainCyclicBoundaryCondition ();
    printHamiltonianAndSpinStatus ();
    cout << "finished to set cyclic boundary condition" << endl;

    //  calculating Hamiltonian
    cout << "start to calculate Hamiltonian" << endl;
    H = 0.0;
    cout << "initialize Hamiltonian" << endl;
    for (i = 1; i < LL - 1; i++) {
        for (j = 1; j < LL - 1; j++) {
            //cout << "hoge" << endl;
            //    cout << "in i = " << i << ", j = " << j << ", s = " << s[i][j] << endl;
            for (k = 0; k < LL; k++) {
                for (l = 0; l < LL; l++) {
                    H += (obtainVanDerWaalsInteractionDiscrete2D (i, j, j, l) +
                        obtainEDLInteractionDiscrete2D (i, j, k, l)) * s[i][j] * s[k][l];
                }
            }
        }
    }
}

double parallelArray2D::obtainProbability (double H) {
    return exp (-1.0 * Beta * H) / Z;

}

void parallelArray2D::obtainCyclicBoundaryCondition () {
    //  cyclic boundary condition 
    for (i = 0; i < LL; i++) s[i][0]      = s[i][LL - 2];
    for (i = 0; i < LL; i++) s[i][LL - 1] = s[i][1];
    for (j = 0; j < LL; j++) s[0][j]      = s[LL - 2][j];
    for (j = 0; j < LL; j++) s[LL - 1][j] = s[1][j];
}


void parallelArray2D::printHamiltonianAndSpinStatus () {
    cout << "Hamiltonian: " << H << endl;
    cout << "Spin status:" << endl;
    for (i = 0; i < LL; i++) {
        for (j = 0; j < LL; j++) {
            cout << s[i][j] << " ";
        }
        cout << endl;
    }
    cout << "End spin status:" << endl;
}



double parallelArray2D::obtainD1DistanceDiscrete1D (int x1, int x2) {
    if (x1 == x2) {
        return 0.0;
    } else {
        return 1.0;
    }
}

double parallelArray2D::obtainD1DistanceDiscrete2D (int x1, int y1, int x2, int y2) {
    if (x1 == x2 && y1 == y2) {
        return 0.0;
    } else {
        return 1.0;
    }
}


double parallelArray2D::obtainEDLInteractionDiscrete2D (int x1, int y1, int x2, int y2) {
    return obtainD1DistanceDiscrete2D (x1, y1, x2, y2) * exp (-1.0 * Kappa * obtainD2DistanceDiscrete2D (x1, y1, x2, y2));

}



double parallelArray2D::obtainVanDerWaalsInteractionDiscrete2D (int x1, int y1, int x2, int y2) {
    return -1.0 * obtainD1DistanceDiscrete2D (x1, y1, x2, y2) * CvdW / (obtainD2DistanceDiscrete2D (x1, y1, x2, y2) + 0.001);

}

double parallelArray2D::obtainD2DistanceDiscrete2D (int x1, int y1, int x2, int y2) {
    return sqrt (pow ((double) (x1 - x2), 2) + pow ((double) (y1 - y2), 2));
}