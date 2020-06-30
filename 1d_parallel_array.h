#include <iostream>
#include <vector>
#include <list>
#include <string>

using namespace std;
class parallelArray1D {
    protected:
        int n_status, n_particles, n_state_of_c;
        int i, j, k, l;
        int c_th;
        double Jedl, JvdW, J;    //  Hamiltonian, interaction terms
        double H, Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        
        vector<int> s, s_record_book, c, c_record_book;
        //vector<string> s_all_state_in_str;
        
        //list<string> all_spin_state;

        void initRandomSpins ();
        void initRandomOrderParameter ();
        void obtainCyclicBoundaryCondition ();
        
        double obtainD1Distance (int i, int j);

        int obtainSpinFromOrderParameter ();

        double obtainEDLInteraction ();
        double obtainVanDerWaalsInteraction ();


        double obtainD1DistanceDiscrete1D (int x1, int x2);    
        double obtainD2DistanceDiscrete1D (int x1, int x2);    
        double obtainEDLInteractionDiscrete1D (int x1, int x2);
        double obtainVanDerWaalsInteractionDiscrete1D (int x1, int x2);
        void executeGibbsSampling ();        
        double obtainConditionalProb1D (int i_particle);
        void transformOrderParameterToSpinState ();
        int obtainSpinFromOrderParameter (int c);
        //double obtainProbability (double H);
    public:
        parallelArray1D ();
                
        void perturbSpinDistributionWithGibbsSampling ();



};