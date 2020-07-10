#include <iostream>
#include <vector>

using namespace std;
class parallelArray2D {
    protected:
        int n_x, n_y, n_paticles, n_spin_state;
        int i_x, i_y, i, j, k, l;
        double D, Dbar, Dtil, Phi, Delta, L, Cexcl;
        double VerPhiTil, VerPhiBar, VerPhi;
        double H, Jedl, JvdW, Alpha;    //  Hamiltonian, interaction terms
        double IonicStrength, Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        
        std::vector <int> s, s_record_book;
        std::vector <double> Sigma;
        //std::vector <vector <int>> s, s_record_book;
        int obtainOneDimPosFromTwoDimPos (int i_x, int i_y); 
        void initRandomSpins ();
        void initSigma ();
        void obtainCyclicBoundaryCondition ();
        double obtainSigmaFromSpin (double s);
        double obtainCondProb1 (int i);
        double obtainCondProb23 (int i);
    public:
        parallelArray2D ();
        
        void printSpinState ();
        
};