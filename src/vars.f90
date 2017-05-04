MODULE globvars

!*************************************************************************************************!
!*
!*        Global Variables Module
!*        
!*  Contains declaration of all global variables and subroutines for:
!*
!*      1) Initialisation of global variables
!*      
!*************************************************************************************************!

  type basisfn                              ! an array of size nbf of this type holds the basis set
    complex(kind=8), dimension(:), allocatable::z        ! Coherent state sub-array - 1 value per dof
    complex(kind=8)                            ::D_big    ! Prefactor - 1 value per basis function
    complex(kind=8), dimension(:), allocatable::d_pes    ! Electronic State pre-exponential - 1 value per state
    real   (kind=8), dimension(:), allocatable::s_pes    ! Classical action - 1 value per electronic state
    complex(kind=8), dimension(:), allocatable::a_pes    ! Electronic state amplitude - 1 value per state and =de^(iS)
  end type basisfn

  type hamiltonian                  ! a matrix of size nbf x nbf of this type hold the hamiltonian
    complex(kind=8), dimension(:,:), allocatable::Hjk      ! matrix of size npes x npes
  end type hamiltonian

  integer:: in_nbf     ! number of basis functions
  integer:: ndim       ! number of dimensions/degrees of freedom
  integer:: npes       ! determines n.o configurations
  integer:: in_pes     ! potential energy surface 1 or 2
  integer:: Ntries     ! Number of recalculations allowed before an error occurs
  integer:: errorflag  ! The global fatal error flag. Stops operation
  integer:: debug      ! Debug flag. Not currently used.
  integer:: conjflg    ! Flag to determine whether Conjugate repetition will occur
  integer:: restrtflg  ! Flag for restarting a previous run for when time has run out
  integer:: reptot     ! Total number of repeats to be carried out in program
  integer:: qsizex     ! size of grid in q,x coordinate
  integer:: psizex     ! size of grid in p,x coordinate
  integer:: qsizey     ! size of grid in q,y coordinate
  integer:: psizey     ! size of grid in p,y coordinate
  integer:: qsizez     ! size of grid in q,z coordinate
  integer:: psizez     ! size of grid in p,z coordinate
  integer:: trainsp    ! number of time steps between adjacent basis functions in a train
  integer:: def_stp    ! default number of basis functions in a train for swarms of trains
  integer:: clonemax   ! maximum number of cloning events allowed per basis function
  integer:: clonefreq  ! minimum number of timesteps between concurrent coloning events

  real(kind=8),external :: ZBQLNOR                ! The normally distributed random number external function
  
  real(kind=8) :: y_00_vp          ! S_0 offset **PSEUDO-TORTION POTENTIAL ONLY**
  real(kind=8) :: A_0_vp           ! S_0 amplitude **PSEUDO-TORTION POTENTIAL ONLY**
  real(kind=8) :: b_0_vp           ! S_0 parabolic parameter **PSEUDO-TORTION POTENTIAL ONLY**
  real(kind=8) :: width_vp         ! S_0 gaussian width **PSEUDO-TORTION POTENTIAL ONLY**
  real(kind=8) :: A_NAC_vp         ! Non-adiabatic Coupling amplitude **PSEUDO-TORTION POTENTIAL ONLY**
  real(kind=8) :: y_10_vp          ! S_1 offset **PSEUDO-TORTION POTENTIAL ONLY**
  real(kind=8) :: b_1_vp           ! S_1 parabolic parameter **PSEUDO-TORTION POTENTIAL ONLY**


  real(kind=8) :: wc_sb            ! Cutoff frequency **SPIN BOSON ONLY**
  real(kind=8) :: beta_sb          ! Temperature parameter **SPIN BOSON ONLY**
  real(kind=8) :: kondo_sb         ! Kondo parameter **SPIN BOSON ONLY**
  real(kind=8) :: wmax_sb          ! Maximum frequency = 5*wc **SPIN BOSON ONLY**
  real(kind=8) :: delta_sb         ! Coupling between PESs in Hamiltonian - provides off diagonals **SPIN BOSON ONLY**
  real(kind=8) :: eps_sb           ! Symmetry parameter used in calculation of diagonals in hamiltonian **SPIN BOSON ONLY**
  
  real(kind=8) :: wc_dl            ! Cutoff frequency **DRUDE LORENTZ ONLY**
  real(kind=8) :: beta_dl          ! Temperature parameter **DRUDE LORENTZ ONLY**
  real(kind=8) :: lambda_dl        ! Reorganisation Energy parameter **DRUDE LORENTZ ONLY**
  real(kind=8) :: wmax_dl          ! Maximum frequency **DRUDE LORENTZ ONLY**
  real(kind=8) :: delta_dl         ! Coupling between PESs in Hamiltonian - provides off diagonals **DRUDE LORENTZ ONLY**
  real(kind=8) :: eps_dl           ! Symmetry parameter used in calculation of diagonals in hamiltonian **DRUDE LORENTZ ONLY**

  real(kind=8) :: mass_fp          ! Mass of a particle **FREE PARTICLE ONLY**

  real(kind=8) :: freq_hp          ! Frequency of a Harmonic Potential **HARMONIC POTENTIAL ONLY**

  real(kind=8) :: dissen_mp        ! Dissasociation Energy **MORSE POTENTIAL ONLY**
  real(kind=8) :: freq_mp          ! frequency of morse potential **MORSE POTENTIAL ONLY**
  real(kind=8) :: mass_mp          ! mass of particle in morse potential **MORSE POTENTIAL ONLY**
  real(kind=8) :: a0_mp            ! Well shape parameter for morse oscillator **MORSE POTENTIAL ONLY**

  real(kind=8) :: mass_iv          ! Particle mass for inverted gaussian **INVERTED GAUSSIAN ONLY**
  real(kind=8) :: freq_iv          ! Laser frequency for inverted gaussian **INVERTED GAUSSIAN ONLY**
  real(kind=8) :: lambda_iv        ! Potential parameter for inverted gaussian **INVERTED GAUSSIAN ONLY**
  real(kind=8) :: inten_iv         ! Laser intensity for inverted gaussian **INVERTED GAUSSIAN ONLY**

  real(kind=8) :: mass_cp          ! Particle mass for Coulomb Potential **COULOMB POTENTIAL ONLY**
  real(kind=8) :: freq_cp          ! Laser frequency for Coulomb Potential **COULOMB POTENTIAL ONLY**
  real(kind=8) :: RC_cp            ! Origin potition for Coulomb Potential **COULOMB POTENTIAL ONLY**
  real(kind=8) :: inten_cp         ! Laser intensity for Coulomb Potential **COULOMB POTENTIAL ONLY**

  real(kind=8) :: lambda_hh        ! Hennon-Heiles Potential Coupling Constant

  real(kind=8) :: sigp             ! Width of the distribution of imaginary parts of the initial wavefunction
  real(kind=8) :: sigq             ! Width of the distribution of real parts of the initial wavefunction
  real(kind=8) :: gam              ! Gamma value, used to calculate the widths sigp and sigq
  real(kind=8) :: initalcmprss     ! Initial compression parameter. Can be automatically tweaked or left static.
  real(kind=8) :: mu               ! Center value for the real and imaginary parts of the initial wavefunction 
  real(kind=8) :: hbar             ! hbar. Currently redundant as hbar assumed to be 1
  real(kind=8) :: Ebfmin           ! Minimum energy allowed for a single basis function
  real(kind=8) :: Ebfmax           ! Maximum energy allowed for a single basis function
  real(kind=8) :: dtmin            ! Minimum stepsize for adaptive step system
  real(kind=8) :: dtmax            ! Maximum stepsize for adaptive step system
  real(kind=8) :: timeend          ! Final time value for simulation
  real(kind=8) :: timestrt         ! Initial time value for simulation (instability occurs when =/= 0)
  real(kind=8) :: dtinit           ! Initial/Only time step size for adaptive/static stepsize
  real(kind=8) :: initsp           ! Grid spacing for regular grids
  real(kind=8) :: bfeps            ! Adaptive basis set cutoff parameter
  real(kind=8) :: uplimnorm        ! upper limit of the norm
  real(kind=8) :: lowlimnorm       ! lower limit of the norm
  real(kind=8) :: sqrtpi           ! square root of pi
  real(kind=8) :: pirl             ! pi

  character(LEN=5) :: ECheck       ! Flag to determine if the energy of basis functions should be checked
  character(LEN=5) :: basis        ! Flag for grids
  character(LEN=5) :: nbfadapt     ! Flag for adaptive basis set size
  character(LEN=5) :: cloneflg     ! Flag for cloning
  character(LEN=5) :: matfun       ! Which matrix calcualtion function (zgesv/zheev) used for linear algebra
  character(LEN=5) :: method       ! MCEv1 or MCEv2 used. Later versions to have other methods available
  character(LEN=2) :: sys          ! System to be simulated. Currently only Spin Boson allowed
  character(LEN=1) :: gen          ! Flag to determine if basis set generation needed (global variant)
  character(LEN=1) :: prop         ! Flag to determine if basis set propagation needed
  character(LEN=1) :: step         ! Flag to determine whether adaptive or static stepsize propagation needed
  character(LEN=1) :: cmprss       ! Flag to determine if automatic tweaking of the compression parameter needed

  complex(kind=8) :: i             ! The imaginary unit

contains

!--------------------------------------------------------------------------------------------------

  subroutine initialise  !  Level 1 Subroutine

    implicit none

    ! This subroutine initialises all global variables to zero and sets the imaginary unit variable.

    in_nbf = 0
    ndim = 0
    npes = 0
    in_pes = 0
    Ntries = 0
    errorflag = 0
    debug = 0
    conjflg = 0
    reptot = 0
    qsizex = 0
    psizex = 0
    qsizey = 0
    psizey = 0
    qsizez = 0
    psizez = 0
    trainsp = 0
    def_stp = 0
    
    y_00_vp = 0.0d0          ! PSEUDO-TORTION POTENTIAL PARAMETER
    A_0_vp = 0.0d0           ! PSEUDO-TORTION POTENTIAL PARAMETER
    b_0_vp = 0.0d0           ! PSEUDO-TORTION POTENTIAL PARAMETER
    width_vp = 0.0d0         ! PSEUDO-TORTION POTENTIAL PARAMETER
    A_NAC_vp = 0.0d0         ! PSEUDO-TORTION POTENTIAL PARAMETER
    y_10_vp = 0.0d0          ! PSEUDO-TORTION POTENTIAL PARAMETER
    b_1_vp = 0.0d0           ! PSEUDO-TORTION POTENTIAL PARAMETER

    wc_sb = 0.0d0            ! SPIN BOSON PARAMETER
    wmax_sb = 0.0d0          ! SPIN BOSON PARAMETER
    kondo_sb = 0.0d0         ! SPIN BOSON PARAMETER
    eps_sb = 0.0d0           ! SPIN BOSON PARAMETER
    delta_sb = 0.0d0         ! SPIN BOSON PARAMETER
    beta_sb = 0.0d0          ! SPIN BOSON PARAMETER
    
    wc_dl = 0.0d0            ! DL SPIN BOSON PARAMETER
    wmax_dl = 0.0d0          ! DL SPIN BOSON PARAMETER
    lambda_dl = 0.0d0        ! DL SPIN BOSON PARAMETER
    eps_dl = 0.0d0           ! DL SPIN BOSON PARAMETER
    delta_dl = 0.0d0         ! DL SPIN BOSON PARAMETER
    beta_dl = 0.0d0          ! DL SPIN BOSON PARAMETER    

    mass_fp = 0.0d0          ! FREE PARTICLE PARAMETER

    freq_hp = 0.0d0          ! HARMONIC OSCILLATOR PARAMETER

    freq_mp = 0.0d0          ! MORSE POTENTIAL PARAMETER
    mass_mp = 0.0d0          ! MORSE POTENTIAL PARAMETER
    dissen_mp = 0.0d0        ! MORSE POTENTIAL PARAMETER
    a0_mp = 0.0d0            ! MORSE POTENTIAL PARAMETER

    mass_iv = 0.0d0          ! INVERTED GAUSSIAN PARAMETER
    freq_iv = 0.0d0          ! INVERTED GAUSSIAN PARAMETER
    lambda_iv = 0.0d0        ! INVERTED GAUSSIAN PARAMETER
    inten_iv = 0.0d0         ! INVERTED GAUSSIAN PARAMETER

    mass_cp = 0.0d0          ! COULOMB POTENTIAL PARAMETER
    freq_cp = 0.0d0          ! COULOMB POTENTIAL PARAMETER
    Rc_cp = 0.0d0            ! COULOMB POTENTIAL PARAMETER
    inten_cp = 0.0d0         ! COULOMB POTENTIAL PARAMETER

    lambda_hh = 0.0d0        ! HENNON-HAILES PARAMETER

    sigp = 0.0d0
    sigq = 0.0d0
    gam = 0.0d0
    initalcmprss = 0.0d0
    mu = 0.0d0
    Ebfmin = 0.0d0
    Ebfmax = 0.0d0
    hbar = 1.0d0
    dtmin = 0.0d0
    dtmax = 0.0d0
    timeend = 0.0d0
    timestrt = 0.0d0
    dtinit = 0.0d0
    initsp = 0.0d0
    bfeps = 0.0d0
    sqrtpi = 1.7724538509055160272981674833411451827975494561223871d0
    pirl = sqrtpi**2.0d0

    i = (0.0d0,1.0d0)

  end subroutine initialise

END MODULE globvars
