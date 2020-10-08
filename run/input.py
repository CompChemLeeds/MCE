#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#					Input File								!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# System (Currently can be SB [Spin Boson], HP [Harmonic Potential], FP [Free Particle], MP [Morse Potential])
#  SB works only with MCE, all others with CCS.
System='SB'

# Spectral Density (only valid for Spin Boson, with values EXP [exponential cutoff], DL [Drude Lorentz/Debye],
#                                or UBO [Underdamped Brownian Oscillator]. LHC [Light-Harvesting complex] to come later)
SpecDen='EXP'

#Restart flag [1/0 for restarting a previously timed out simulation. This should be automatically set by the restart script]
####restart=0

#Flag for adaptive altering of the compression parameter in swarms/grid-swarms/train-swarms or the grid spacing for grids. [YES/NO]
cmprss='YES'

#Propagation method [can be CCS, MCEv1, MCEv2, or MCE12 (uses both MCE methods)]
method='MCEv2'

#Conjugate repeats flag. Allows even number repeats to start at the complex conjugate of the previous runs initial position [YES/NO]
Conjugate_Repeats='YES'

# Number of basis functions
in_nbf=10

# Number of dimensions
ndim=20

#determine which PES you are on initially
in_PES=1

#determine n.o of configurations for MCE
npes=2

# determines shape of the initial basis [SWARM/SWTRN]
basis='SWARM'

# allows for the use of the quantum superposition sampling amplitudes [1/0]. 0=standard initial sampling, 1=quantum superposition sampling.
qss=1

# Compression parameter (can be tweaked to allow better norm value, or altered automatically)   Used for swarms/swarm-trains
ALCMP=1.0

# frequency flag for using a set of precalculated frequencies with the spin boson model
# freqflg = 1 enables reading of precalculated frequencies, freqflg = 0 mean they are calculated withing the program
freqflg=0

# Flag for cloning basis functions (yes/no/blind/blind+)
Cloning='no'

# Cloning threshold (value of |sum_r(a_{r,k})|) - must be >= 0.05 and < 0.25, default 0.249
Threshold=0.249

# Maximum number of Cloning events allowed
max_cloning=4

# Minimum cloning frequency (ie how many timesteps since last cloning is new cloning event allowed)
clon_freq=200

#Quantum Superposition Cloning exclusion paramter between the two child trajectories should >= ??? and < ???
QSC_epsilon=0.1 

# spacing between trajectories for train type basis set. Used for swarm-trains
trainsp=300

# length of train in carriages. Used for swarm-trains only
train_len=10

# Size of the central swarm for swtrn basis. Used for swarm-trains only
swtrn_swarm_size=10

# Random Number generation function (ZBQL - using ZBQLNOR subroutine, GAUS - using function based on numerical recipes)
randfunc='ZBQL'

# Seed value for doing the random number routine- if you do not specify it (leave default value of 0) will automatically generate one
SEED=0

# gamma factor
gamma=1.0

# center of initial random gaussian
mu=0.0

# hbar if left commented or blank, defaults to 1.0d0
hbar=1.0

# end of input
