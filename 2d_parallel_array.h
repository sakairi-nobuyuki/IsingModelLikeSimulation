#include <iostream>
#include <vector>

using namespace std;
class parallelArray2D {
    protected:
        int n_x, n_y, n_particles, n_spin_state;
        int s_gen[3][2] = {{1, 0}, {0, 1}, {1, 1}}, e_x[2] = {1, 0}, e_y[2] = {0, 1};    //  generator of spin with 2 degree of freedom.
        int n_particles_perturb;
        int i_x, i_y, i, j, k, l, i_sample;
        double D, Dbar, Dtil, Phi, Delta, L, Cexcl;
        double VerPhiTil, VerPhiBar, VerPhi;
        double H, Htmp, Jedl, JvdW, Alpha;    //  Hamiltonian, interaction terms
        double IonicStrength, Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        double Tsa, TsaMax, TsaMin, dTsa, ksa, Psa, PsaRef;
        
        std::vector <int> s, s_record_book, s_tmp;
        std::vector <double> VerPhiEySij, VerPhiEySneib, VerPhiExSij, VerPhiExSneib;
        std::vector <double> Sigma;
        //std::vector <vector <int>> s, s_record_book;
        int obtainOneDimPosFromTwoDimPos (int i_x, int i_y); 
        void initRandomSpins ();
        void obtainCyclicBoundaryCondition ();

        
        void applyCyclicBoundaryCondition (vector<int>& spin);        
        double obtainHamiltonian (vector<int>& spin);
        //double obtainHamiltonian (vector<int>& spin);
    public:
        parallelArray2D ();
        
        void executeSimulatedAnnealing ();
        void printSpinState ();
        void obtainJointProbability ();
};