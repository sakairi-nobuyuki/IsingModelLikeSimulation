#  number of particles in x-direction
n_x = [20]            
#  number of particles in y-direction
n_y = [20]            
#  degree of freedom of the spin
n_spin_state = [3]    

###   physical chemicstry parameters
#  salt concentration (mol/L)
IonicStrength = [0.001, 0.01, 0.1]   
#  volume fraction  (-)
Phi   = [0.01, 0.05, 0.1]     
#  width of the particle (nm)
Delta = [1.0e-09]     
#  length of  th particle (nm)
L     = [500.0e-09]           
#  coefficient of exclusive volume (-)
Cexcl = [0.8]                 

###   Thermodynamic parameters
#  Boltzmann factor (J/K)
kB    = [1.0]                 
#  Temperature (K) 
T     = [100, 300, 500]       

###   Simulated Annealing parameters    
#  max Temperature    
TsaMax = [1.0]
#  minimum temperature
TsaMin = [0.1]
#  temperature step
dTsa   = [0.001]
#  Apparent Boltzmann factor
ksa    = [1.0e-23]   
#  reference probability
PsaRef = [0.1]    
#  trial at one annealing stage
n_trial  = [10]
###  initialize geometrical parameters
# #  ratio of van der Waals force against EDL force (no meaning right now)
CvdW     = [0.0001]
#  EDL force coefficient (no meaning right now)
Alpha    = [1.0e+09]

###  Output setting
# output file name base
output_file_name = ['out']   
# output frequency in every annealing stage
output_freq      = [100]              