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
        double Phi, Delta, L, D, Dex, Cex, IonicStrength;   //  Dex: distance concerning to exclusive volume, Cex: coefficient concerning to exclusive voume
        double Alpha;   //  Ratio of EDL force and van der Waals force
        double H, Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double ProtoProb;
        double Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        
        vector<int> s, s_record_book, c, c_record_book;

        void initRandomOrderParameter ();
        void printSpinAndOrderParameterState ();
        void obtainCyclicBoundaryCondition ();
        
        double obtainD1Distance (int i, int j);

        int obtainSpinFromOrderParameter ();

        double obtainEDLInteraction ();
        double obtainVanDerWaalsInteraction ();


        double obtainD1DistanceDiscrete1D (int x1, int x2);    
        double obtainD2DistanceDiscrete1D (int x1, int x2);    
        double obtainEDLInteractionDiscrete1D (int x1, int x2);
        double obtainVanDerWaalsInteractionDiscrete1D (int x1, int x2);
  
        double obtainInteractionPotential (int s_1, int s_2);
        double obtainConditionalProb1D (int i_particle, int c_subject);
        void transformOrderParameterToSpinState ();
        int obtainSpinFromOrderParameter (int c);
        
    public:
        parallelArray1D ();

        void executeGibbsSampling ();              
        void perturbSpinDistributionWithGibbsSampling ();



};