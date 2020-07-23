#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <string>
#include <numeric>

#include "2d_parallel_array.h"

using namespace std;

parallelArray2D::parallelArray2D () {
    loadConfigFile ();
    cout << "load config file" << endl;
    //n_spin_state = 3;
    //n_x   = 10;
    //n_y   = 10;
    n_particles = n_x * n_y;
    
    //  allocate spin and spin record book
    s.resize (n_particles);
    s_tmp.resize (n_particles);
    s_record_book.resize (n_particles);
    Sigma.resize (n_particles);
    VerPhiEySij.resize (n_particles);
    VerPhiEySneib.resize (n_particles);
    VerPhiExSij.resize (n_particles);
    VerPhiExSneib.resize (n_particles);

    //IonicStrength = 0.01;
    //Phi   = 0.05;        //  volume fraction
    //Delta = 1.0e-09;     //  width of the particle
    //L     = 500.0e-09;   //  length of  th particle
    //Cexcl = 0.8;

    //kB = 1.0;
    //T  = 100;
    Beta = 1.0 / (kB * T);
    
    //TsaMax = 1.0;
    //TsaMin = 0.1;
    //dTsa   = 0.001;

    //  initialize geometrical parameters
    //Kappa = 1.0 / (0.304 / sqrt (IonicStrength) * 1.0e-09);   
    Kappa = 1.0 / (0.304 / sqrt (IonicStrength) * 1.0e-09);   
    //Dbar = Delta * (1.0 / Phi - 1.0);    ///  1 dimensional case
    Dbar = sqrt (2.0 * L * Delta / Phi) - 0.5 * Delta;
    Dtil = 0.5 * (L - Delta) * Delta / L * Cexcl;
    //CvdW  = 0.0001;
    //Alpha = 1.0e+09;
    
    VerPhiBar = exp (-Kappa * Dbar);
    VerPhiTil = 1.0 - exp (Kappa * Dtil);
    VerPhi    = VerPhiBar * (pow (VerPhiTil, 2) - VerPhiTil);

    //ksa = 1.0e-23;
    //PsaRef = 0.1;    

    cout << "initializin geometrical parameters" << endl;
    cout << "  phi:   " << Phi << "(-)" << endl;
    cout << "  kappa: " << Kappa * 1.0e-09 << "(1/m-9)" << endl;
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
    cout << "  s_u1:    " << "[" << s_gen[0][0] << ", " << s_gen[0][1] << "]" << endl;
    cout << "  s_u2:    " << "[" << s_gen[1][0] << ", " << s_gen[1][1] << "]" << endl;
    cout << "  s_u3:    " << "[" << s_gen[2][0] << ", " << s_gen[2][1] << "]" << endl;

    // initialization with random spins
    initRandomSpins ();
    applyCyclicBoundaryCondition (s);
    H = obtainHamiltonian (s);
    //initSigma ();

//    loadConfigFile ();
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




void  parallelArray2D::executeSimulatedAnnealing () {
    int n_shifted_particles;
    for (Tsa = TsaMax; Tsa > TsaMin; Tsa -= dTsa) {
        //  obtain number of particles to make perturbation
        n_particles_perturb = (int) (n_particles * (Tsa - TsaMin) / (TsaMax - TsaMin));
        cout << "n_particles to perturb: " << n_particles_perturb << endl;
        //  choose a set of particles to perturbate
        for (i = 0; i < n_particles; i++) s_tmp[i] = s[i];
        for (i = 0; i < n_particles; i++) s_record_book[i] = 0;
        for (i = 0; n_shifted_particles = accumulate (s_record_book.begin (), s_record_book.end (), 0) < n_particles_perturb; i++) {
            i_sample = rand () % n_particles;
            if (s_record_book[i_sample] == 0) {
                s_tmp[i_sample] = rand () % n_spin_state;
                s_record_book[i_sample] = 1;
            }
        }

        //  obtain new Hamiltonian
        //for (i = 0; i < n_particles; i++) {
        //    if (s_tmp[i] != s[i]) cout << "i: " << i << " s_tmp: " << s_tmp[i] << ", s: " << s[i] << endl;
        //}
        applyCyclicBoundaryCondition (s_tmp);
        Htmp = obtainHamiltonian (s_tmp);

        //  evaluate which cases shall be applied
        cout << "H = " << H << " Htmp = " << Htmp << " H - Htmp = " << H - Htmp << endl;
        
        Psa = exp ((H - Htmp) / (ksa * Tsa));
        cout << "Psa = " << Psa << " PsaRef = " << PsaRef << endl;
        if (Htmp < H) {
            for (i = 0; i < n_particles; i++) s[i] = s_tmp[i];
            H = Htmp;
        } else if (Psa > PsaRef) {
            for (i = 0; i < n_particles; i++) s[i] = s_tmp[i];
            H = Htmp;
        }
        



    }
}


void parallelArray2D::obtainJointProbability () {
    H = obtainHamiltonian (s);
    cout << "probability: " << exp (-Beta * H);

}

double parallelArray2D::obtainHamiltonian (vector<int>& spin) {
    double Hamiltonian;
    //for (i = 0; i < n_paticles; i++) {
    cout << "going to obtain Hamiltonian" << endl;
    applyCyclicBoundaryCondition (spin);
    cout << "finished cyclic Boundary Condition" << endl;
    Hamiltonian = 0.0;
    //for (i = n_x; i < n_particles - n_x; i++) {
    for (j = 1; j < n_y - 1; j++) {
        for (i = 1; i < n_x - 1; i++) {
            Hamiltonian += 
                (1.0 - VerPhiTil * s_gen[spin[obtainOneDimPosFromTwoDimPos (i, j)]][1]) 
                    * (2.0 - VerPhiTil * (s_gen[spin[obtainOneDimPosFromTwoDimPos (i - 1, j)]][1]) 
                    *                   + s_gen[spin[obtainOneDimPosFromTwoDimPos (i + 1, j)]][1])
                + (1.0 - VerPhiTil * s_gen[spin[obtainOneDimPosFromTwoDimPos (i, j)]][0]) 
                    * (2.0 - VerPhiTil * (s_gen[spin[obtainOneDimPosFromTwoDimPos (i, j - 1)]][0]) 
                    *                   + s_gen[spin[obtainOneDimPosFromTwoDimPos (i, j + 1)]][0]);
            //cout << "i: " << obtainOneDimPosFromTwoDimPos (i, j);
            //cout << ", i neib " <<  obtainOneDimPosFromTwoDimPos (i-1, j) << ", " << obtainOneDimPosFromTwoDimPos (i+1, j);
            //cout << ", " << obtainOneDimPosFromTwoDimPos (i, j-1) << ", " << obtainOneDimPosFromTwoDimPos (i, j+1) << endl;
            //cout << " s: " << spin[obtainOneDimPosFromTwoDimPos (i, j)];
            //cout << ", s neib: " << spin[obtainOneDimPosFromTwoDimPos (i-1, j)] << ", " << spin[obtainOneDimPosFromTwoDimPos (i+1, j)];
            //cout << ", " << spin[obtainOneDimPosFromTwoDimPos (i, j-1)] << ", " << spin[obtainOneDimPosFromTwoDimPos (i, j+1)]<< endl;
            //cout << " e_y s i-1 j: " << s_gen[spin[obtainOneDimPosFromTwoDimPos (i-1, j)]][1];
            //cout << ", e_y s i+1 j: " << s_gen[spin[obtainOneDimPosFromTwoDimPos (i+1, j)]][1] << endl;
            //cout << " e_x s i j-1: " << s_gen[spin[obtainOneDimPosFromTwoDimPos (i, j-1)]][0];
            //cout << ", e_x s i j-1: " << s_gen[spin[obtainOneDimPosFromTwoDimPos (i, j+1)]][0] << endl;
            //cout << " H: " << Hamiltonian << endl;
        }
    }
    Hamiltonian *= VerPhiBar;
    cout << "Hamiltonian: " << Hamiltonian << endl;

    return Hamiltonian;
}

void parallelArray2D::applyCyclicBoundaryCondition (vector<int>& spin) {
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
    //printSpinState ();
}



inline int parallelArray2D::obtainOneDimPosFromTwoDimPos (int i_x, int i_y) {
    return i_y * n_x + i_x;

}


void parallelArray2D::loadConfigFile () {
    fstream filestream ("conf_2d.dat");
    string buff, token, comment_reduced_buff, separated_buff;
    
    int i_token, i_lines;

    //if (!filestream.is_open ())   return 0;

    for (i_lines = 0; !filestream.eof (); i_lines++) {
        getline (filestream, buff);
        if (buff[0] == '#')  continue;
        
        istringstream comment_reduced_stream (buff);
        getline (comment_reduced_stream, comment_reduced_buff, '#');
        //cout << buff << endl;
        //cout << comment_reduced_buff << endl;

        istringstream separated_stream (comment_reduced_buff);
        vector<string> separated_string;    
        for (i_token = 0; getline (separated_stream, separated_buff, ':'); i_token++) {
            separated_string.push_back (separated_buff);
        }

        if (separated_string[0] == "n_x") {
            //cout << separated_string[0] << "=" << separated_string[1] << endl;
            //cout << typeid (separated_string[1]).name () << endl;
            //n_x = static_cast<int> (separated_buff[1]);
            //n_x = std::atoi (separated_string[1].c_str);
            n_x = stoi (separated_string[1]);
            //cout << n_x << " cast " << std::stoi ("20") << endl;
        }
        if (separated_string[0] == "n_y") {
            n_y = stoi (separated_string[1]);
            //cout << n_y << endl;
        }
        if (separated_string[0] == "n_spin_state") {
            n_spin_state = stoi (separated_string[1]);
            //cout << n_spin_state << endl;
        }
        if (separated_string[0] == "IonicStrength") {
            IonicStrength = stod (separated_string[1]);
            //cout << separated_string[1] << endl;
            //cout << IonicStrength << endl;
        }
        if (separated_string[0] == "Phi")   Phi   = stod (separated_string[1]);
        if (separated_string[0] == "Delta") {
            Delta = stod (separated_string[1]);
            //cout << Delta << endl;
        }
        if (separated_string[0] == "L")     L     = stod (separated_string[1]);
        if (separated_string[0] == "Cexcl") Cexcl = stod (separated_string[1]);
        if (separated_string[0] == "kB")    kB    = stod (separated_string[1]);
        if (separated_string[0] == "T")     T    = stod (separated_string[1]);
        if (separated_string[0] == "TsaMax")     TsaMax    = stod (separated_string[1]);
        if (separated_string[0] == "TsaMin")     TsaMin    = stod (separated_string[1]);
        if (separated_string[0] == "dTsa")       dTsa      = stod (separated_string[1]);
        if (separated_string[0] == "ksa")        ksa       = stod (separated_string[1]);
        if (separated_string[0] == "PsaRef")     PsaRef    = stod (separated_string[1]);
        if (separated_string[0] == "CvdW")       CvdW      = stod (separated_string[1]);
        if (separated_string[0] == "Alpha") {
            Alpha     = stod (separated_string[1]);
            //cout << Alpha << endl;
        }   
        
    }
    
    
}