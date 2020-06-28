#include <iostream>
#include <vector>
#include <list>
#include <string>

using namespace std;
class parallelArray1D {
    protected:
        int n_status, n_particles;
        int i, j, k, l;
        double Jedl, JvdW, J;    //  Hamiltonian, interaction terms
        double H, Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double Z, Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        
        vector<int> s, s_record_book;
        vector<string> s_all_state_in_str;
        
        //list<string> all_spin_state;

        void initRandomSpins ();
        void obtainCyclicBoundaryCondition ();
        
        double obtainD1Distance (int i, int j);


        void setAllSpinState (list<string>& AllSpinState);

        double obtainEDLInteraction ();
        double obtainVanDerWaalsInteraction ();


        double obtainD1DistanceDiscrete1D (int x1, int x2);    
        double obtainD2DistanceDiscrete1D (int x1, int x2);    
        double obtainEDLInteractionDiscrete1D (int x1, int x2);
        double obtainVanDerWaalsInteractionDiscrete1D (int x1, int x2);
        double obtainDistributionFunction ();
        double obtainProbability (double H);
    public:
        parallelArray1D ();
        
        //double obtainHamiltonian (int *s, int n_particles);
        double obtainHamiltonian (vector<int>& s);
        void printHamiltonianAndSpinStatus ();
        void perturbSpinDistributionWithGibbsSampling ();



};