
!--------------------------------------------------------------------------

MODULE intvars

  integer:: ndim       ! number of dimensions/degrees of freedom
  integer :: freqflg            ! Flag for reading pre-calculated frequencies (1=read, 0=calculate)
  integer :: errorflag
    
  real(kind=8) :: wc_exp            ! Cutoff frequency         *Exp Cutoff*
  real(kind=8) :: kondo_exp         ! Kondo parameter          *Exp Cutoff*
  real(kind=8) :: wmax_exp          ! Maximum frequency = 5*wc *Exp Cutoff*
  
  real(kind=8) :: wc_dl            ! Cutoff frequency                **DRUDE LORENTZ ONLY**
  real(kind=8) :: lambda_dl        ! Reorganisation Energy parameter **DRUDE LORENTZ ONLY**
  real(kind=8) :: wmax_dl          ! Maximum frequency               **DRUDE LORENTZ ONLY**
  
  real(kind=8) :: lambda_ubo       ! Reorganisation Energy **Underdamped Brownian Oscillator**
  real(kind=8) :: gamma_ubo        ! Peak width            **Underdamped Brownian Oscillator**
  real(kind=8) :: w0_ubo           ! Peak Frequency        **Underdamped Brownian Oscillator**
  real(kind=8) :: wmax_ubo         ! Maximum Frequency     **Underdamped Brownian Oscillator**

  real(kind=8) :: s                ! exponent for spectral densities - s<1 = sub-ohmic, s=1 = ohmic, s>1 = super-ohmic
  real(kind=8) :: pi_rl            ! pi
  
  character(LEN=3) :: specden      ! the spectral density type (EXP,DL or UBO ---- LHC to come)

contains

!--------------------------------------------------------------------------------------------------

  subroutine int_initialise  !  Level 1 Subroutine

    implicit none

    ! This subroutine initialises all global variables to zero and sets pi

    ndim = 0
    freqflg = 0
    errorflag = 0
    
    wc_exp = 0.0d0           ! Cutoff frequency   *Exp Cutoff*
    kondo_exp = 0.0d0        ! Kondo parameter    *Exp Cutoff*
    wmax_exp = 0.0d0         ! Maximum frequency  *Exp Cutoff*
  
    wc_dl = 0.0d0            ! Cutoff frequency                 **DRUDE LORENTZ ONLY**
    lambda_dl = 0.0d0        ! Reorganisation Energy parameter  **DRUDE LORENTZ ONLY**
    wmax_dl = 0.0d0          ! Maximum frequency                **DRUDE LORENTZ ONLY**
  
    lambda_ubo = 0.0d0       ! Reorganisation Energy **Underdamped Brownian Oscillator**
    gamma_ubo = 0.0d0        ! Peak width            **Underdamped Brownian Oscillator**
    w0_ubo = 0.0d0           ! Peak Frequency        **Underdamped Brownian Oscillator**
    wmax_ubo = 0.0d0         ! Maximum Frequency     **Underdamped Brownian Oscillator**

    s = 0.0d0                ! exponent for spectral densities - s<1 = sub-ohmic, s=1 = ohmic, s>1 = super-ohmic
    pi_rl = 3.1415926535897932384d0
    
    write(0,*) "All parameters initialised"

  end subroutine int_initialise

END MODULE intvars

!************************************************************************************************************************!

MODULE calculations

  use intvars

contains 

  subroutine spectral_density (wm,Jw,j,an)
  
    implicit none
    
    real(kind=8) :: wc, lambda_temp
    real(kind=8), intent(in) :: wm
    real(kind=8), intent(inout) :: Jw
    integer, intent(in) :: j, an
    
    if (errorflag .ne. 0) return

    select case (specden)
    
      case("EXP")
         
        lambda_temp = pi_rl*kondo_exp*wc_exp/(4.0d0*s)
    
        if (an .eq. 1) then
          Jw= 2.0d0*lambda_temp*((wm/wc_exp)**s)
        else
          Jw= 2.0d0*s*lambda_temp*(wm**(s-1.0d0))*(wc_exp**(-1.0d0*s))
          Jw = Jw * dexp(-1.0d0*wm/wc_exp)
        end if
        
      case("DL")
            
        if (an .eq. 1) then
          Jw= 2.0d0*lambda_dl*((wm/wc_dl)**s)
        else
          Jw= 2.0d0*s*lambda_dl*(wm**(s-1.0d0))*(wc_dl**(-1.0d0*s))
          Jw = Jw * (1.0d0/(1.0d0+((wm/wc_dl)**2.0d0)))
        end if
        
      case("UBO")
          
        wc = (w0_ubo**2.0d0)/gamma_ubo 
 
        if (an .eq. 1) then
          Jw= 2.0d0*lambda_ubo*((wm/wc)**s)
        else
          Jw= 2.0d0*s*lambda_ubo*(wm**(s-1.0d0))*(wc**(-1.0d0*s))
          Jw = Jw * (w0_ubo**4.0d0)/((((w0_ubo**2.0d0)-(wm**2.0d0))**2.0d0)+((gamma_ubo**2.0d0)*(wm**2.0d0)))
        end if
        
      case default
      
        write(0,*) "Error! The spectral density is not recognised! Aborting"
        write(0,*) "Got a value of ", specden
        errorflag = 1
      
    end select
    
    return
    
  end subroutine spectral_density
  
!--------------------------------------------------------------------------

  subroutine lhc (wm,Jw,s,j,an)
  
    implicit none
    
    real(kind=8) :: lambda, w01, w02, w03, gamm1, gamm2, gamm3, wc, j1, j2, j3, j4
    real(kind=8) :: w04, w05, w06, w07, gamm4, gamm5, gamm6, gamm7, j5, j6, j7, j8
    real(kind=8), intent(in) :: wm, s
    real(kind=8), intent(inout) :: Jw
    integer, intent(in) :: j, an
    
    if (errorflag .ne. 0) return
    
    lambda = 2.0d0
    wc = 5.0d0
    w01 = 10.0d0
    w02 = 15.0d0
    w03 = 20.0d0
    w04 = 25.0d0
    w05 = 30.0d0
    w06 = 35.0d0
    w07 = 40.0d0
    gamm1 = 0.25d0
    gamm2 = 0.25d0
    gamm3 = 0.25d0
    gamm4 = 0.25d0
    gamm5 = 0.25d0
    gamm6 = 0.25d0
    gamm7 = 0.25d0
    
    
    if (s.ne.1.0d0) then
      write(6,*) "Can only do LHC with s=1"
      stop
    end if
    
    if (an .eq. 1) then
      Jw= 2.0d0*lambda*(wm/wc)
    else
      j1 = 2.0d0*lambda*wc/((wm**2.0d0)+(wc**2.0d0))
      j2 = gamm1**2.0d0*w01/(((w01**2.0d0-wm**2.0d0)**2.0d0)+(gamm1**2.0d0*wm**2.0d0))
      j3 = gamm2**2.0d0*w02/(((w02**2.0d0-wm**2.0d0)**2.0d0)+(gamm2**2.0d0*wm**2.0d0))
      j4 = gamm3**2.0d0*w03/(((w03**2.0d0-wm**2.0d0)**2.0d0)+(gamm3**2.0d0*wm**2.0d0))
      j5 = gamm4**2.0d0*w04/(((w04**2.0d0-wm**2.0d0)**2.0d0)+(gamm4**2.0d0*wm**2.0d0))
      j6 = gamm5**2.0d0*w05/(((w05**2.0d0-wm**2.0d0)**2.0d0)+(gamm5**2.0d0*wm**2.0d0))
      j7 = gamm6**2.0d0*w06/(((w06**2.0d0-wm**2.0d0)**2.0d0)+(gamm6**2.0d0*wm**2.0d0))
      j8 = gamm7**2.0d0*w07/(((w07**2.0d0-wm**2.0d0)**2.0d0)+(gamm7**2.0d0*wm**2.0d0))
      Jw = j1 + j2 + j3 + j4 + j5 + j6 + j7 + j8 
    end if
    
    return
    
  end subroutine lhc
  
!--------------------------------------------------------------------------------------------------

  subroutine readspecparams  
  
    implicit none
    character(LEN=100)::LINE
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=10028, file='inham.dat', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening inham.dat file'
      errorflag = 1
      return
    end if

    read(10028,*,iostat=ierr)LINE

    do while (ierr==0)

      if (LINE=='SB_s_factor') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,s
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading 's' exponent value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBEwc') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,wc_exp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wc value for exponential cutoff"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBEkondo') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,kondo_exp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading kondo parameter value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBEwmaxfact') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,wmax_exp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wmax value for exponential cutoff"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBDwc') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,wc_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wc value for drude lorentz cutoff"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBDlambda') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,lambda_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading lambda value for drude lorentz cutoff"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBDwmaxfact') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,wmax_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wmax factor value for drude lorentz cutoff"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBBw0') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,w0_ubo
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading w0 value for underdamped brownian oscillator"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBBlambda') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,lambda_ubo
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading lambda value for underdamped brownian oscillator"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBBGamma') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,gamma_ubo
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading gamma parameter value for underdamped brownian oscillator"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='SBBwmaxfact') then
        backspace(10028)
        read(10028,*,iostat=ierr)LINE,wmax_ubo
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wmax factor value for underdamped brownian oscillator"
          errorflag = 1
          return
        end if
        n = n+1
      end if

      read(10028,*,iostat=ierr) LINE

    end do
    
    close (10028)
    
    if (n.ne.11) then
      write(0,"(a,i0)") "Not all required variables read in readspecparams subroutine. n = ", n
      errorflag = 1
      return
    end if
  
  end subroutine readspecparams
  
!----------------------------------------------------------------------------------------------

  subroutine readintparams  
  
    implicit none
    character(LEN=100)::LINE
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=10027, file='input.dat', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening input.dat file'
      errorflag = 1
      return
    end if

    read(10027,*,iostat=ierr)LINE

    do while (ierr==0)

      if (LINE=='SpecDen') then
        backspace(10027)
        read(10027,*,iostat=ierr)LINE,specden
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading spectral density value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='ndim') then
        backspace(10027)
        read(10027,*,iostat=ierr)LINE,ndim
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading ndim value"
          errorflag = 1
          return
        end if
        n = n+1  
      else if (LINE=='freqflg') then
        backspace(10027)
        read(10027,*,iostat=ierr)LINE,freqflg
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading freqflg value"
          errorflag = 1
          return
        end if
        n = n+1
      end if
      
      read(10027,*,iostat=ierr)LINE
      
    end do
    
    close (10027)
      
    if (ndim.le.0) then
      write(0,"(a)") "Number of Degrees of Freedom <= 0"
      ierr=-1
      errorflag = 1
      return
    end if
    
    if (n.ne.3) then
      write(0,"(a,i0)") "Not all required variables read in readintparams subroutine. n = ", n
      errorflag = 1
      return
    end if
    
  end subroutine readintparams
  
end module calculations  

!************************************************************************************************************

program integrator

  use intvars
  use calculations

  implicit none
  real(kind=8), dimension (:), allocatable :: wm, a, x, Cm, summation, freq
  real(kind=8) :: wmax, wtest, del_w, fact, tar, Jw
  integer :: i, j, m, points, sumpoints, ierr, low, high, mid, an, block
  
  call int_initialise
  call readspecparams
  call readintparams
  
  if (specden=="EXP") then 
    wmax = wc_exp * wmax_exp
  else if (specden=="DL") then 
    wmax = wc_dl * wmax_dl
  else if (specden=="UBO") then 
    wmax = w0_ubo * wmax_ubo
  else
    errorflag = 1
  end if
  
  points = 500000
  sumpoints=points/2
  
  block=5000
  del_w = wmax/points
  
  if (s.ge.1.0) then
    block = 0
  end if
  
  allocate (a(0:points))
  allocate (summation(1:sumpoints))
  allocate (freq(1:sumpoints))
  
  summation = 0.0d0
  a = 0.0d0
  an = 0
  
  do j=0,points
    if (errorflag .ne. 0) exit
    wtest = dble(j)*del_w
    call spectral_density(wtest, a(j), j, an)
!    call lhc(wtest, a(j), j, an)
  end do
  
  do j=1,sumpoints
    if (errorflag .ne. 0) exit
    if (j.le.block) then
      wtest = dble(j)*del_w*2.0d0
      i=j*2
      an=1
      call spectral_density(wtest, summation(j), i, an)
    else
      summation(j) = summation(block) + a(block*2) +a(2*j)
      do i=1,((2*j)-1)
        summation(j) = summation(j) + 2.0d0*a(i)
        if (mod(i,2) .eq. 1) then
          summation(j) = summation(j) + 2.0d0*a(i)
        end if
      end do
    end if
    summation(j) = summation(j) * (del_w/3.0d0)
    freq(j) = dble(j)*del_w*2.0d0
    an=0
  end do
      
  fact = summation(sumpoints)/dble(ndim)
    
  allocate (wm(ndim))
  allocate (Cm(ndim))
   
  do m=1,ndim
    if (errorflag .ne. 0) exit
    tar = m*fact
    low = 0
    high = sumpoints+1
    
    do while (high-low.gt.1)
      if (errorflag .ne. 0) exit
      mid = (high+low)/2
      if ((summation(sumpoints).ge.summation(1)).eqv.(tar.ge.summation(mid))) then
        low = mid
      else
        high = mid
      end if
    end do
    if (tar.eq.summation(1)) then
      j = 1
    else if (tar.eq.summation(sumpoints)) then
      j = sumpoints - 1
    else
      j = low
    end if
    
    j=j+1
    
    wm(m) = freq(j)

  end do
  
  if (errorflag .eq. 0) then
    
    fact = sqrt(fact * 2.0d0/pi_rl)
  
    do m=1,ndim
      Cm(m) = wm(m)*fact
    end do
  
    open(unit=100, file="freq.dat", status="unknown", iostat=ierr)
    
    if (ierr.eq.0) then
    
      do m=1,ndim
        write (100,"(i0,2(1x,es16.8e3))")  m, wm(m), Cm(m)
      end do
  
      close (100)
      
    end if
    
  end if
    
end program integrator

