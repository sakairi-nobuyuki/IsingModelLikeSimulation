#include <iostream>
#include <vector>

using namespace std;
class parallelArray2D {
    protected:
        int LL;
        int i, j, k, l;
        double H, Jedl, JvdW;    //  Hamiltonian, interaction terms
        double Kappa, CvdW, d;   //  width of EDL, coeffient of van der Waals force
        double Z, Beta, kB, T;      //  Beta = kB T, kB: Boltzmann factor, T: temperature
        
        std::vector <vector <int>> s, s_record_book;

        void initRandomSpins ();
        void obtainCyclicBoundaryCondition ();
        void perturbSpinDistributionWithGibbsSampling ();
        double obtainD1Distance (int i, int j);

        double obtainEDLInteraction ();
        double obtainVanDerWaalsInteraction ();


        double obtainD1DistanceDiscrete1D (int x1, int x2);    
        double obtainD1DistanceDiscrete2D (int x1, int y1, int x2, int y2);
        double obtainEDLInteractionDiscrete2D (int x1, int y1, int x2, int y2);
        double obtainVanDerWaalsInteractionDiscrete2D (int x1, int y1, int x2, int y2);
        double obtainD2DistanceDiscrete2D (int x1, int y1, int x2, int y2);
        double obtainDistributionFunction ();
        double obtainProbability (double H);
    public:
        parallelArray2D ();
        
        void obtainHamiltonian ();
        void printHamiltonianAndSpinStatus ();



};