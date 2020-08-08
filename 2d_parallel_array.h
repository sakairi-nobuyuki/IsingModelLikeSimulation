#include <iostream>
#include <vector>
#include <random>

using namespace std;
class parallelArray2D {
    protected:
        int n_x, n_y, n_particles, n_spin_state;
        int s_gen[3][2] = {{1, 0}, {0, 1}, {1, 1}}, e_x[2] = {1, 0}, e_y[2] = {0, 1};    //  generator of spin with 2 degree of freedom.
        int n_particles_perturb;
        int i_x, i_y, i, j, k, l, i_sample, i_trial, n_trial, n_tsa;
        int output_freq;
        double D, Dbar, Dtil, Phi, Delta, L, Cexcl;
        double VerPhiTil, VerPhiBar, VerPhi;
        double H, Htmp, Jedl, JvdW, Alpha;    //  Hamiltonian, interaction terms
        double IonicStrength, Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        double Tsa, TsaMax, TsaMin, dTsa, ksa, Psa, PsaRef;
        string output_file_name;
        
        std::vector <int> s, s_record_book, s_tmp;
        std::vector <double> VerPhiEySij, VerPhiEySneib, VerPhiExSij, VerPhiExSneib;
        std::vector <double> Sigma;

        std::random_device prd;
        mt19937_64 mt_64;
        uniform_int_distribution<int> prd_spin_state, prd_n_particles;

        //std::vector <vector <int>> s, s_record_book;

        int loadConfigFile ();

        int obtainOneDimPosFromTwoDimPos (int i_x, int i_y); 
        void initRandomSpins ();
        void initKsa ();
        void initPseudoRandomNumber ();
        void obtainCyclicBoundaryCondition ();

        
        void applyCyclicBoundaryCondition (vector<int>& spin);        
        double obtainHamiltonian (vector<int>& spin);
        //double obtainHamiltonian (vector<int>& spin);
        void renewHamiltonianAndSpinState ();
        int initOutputFile ();
        int dumpSpinState (int n);
    public:
        parallelArray2D ();
        
        void executeSimulatedAnnealing ();
        void printSpinState ();
        void obtainJointProbability ();

        int writeOutputFile ();
};