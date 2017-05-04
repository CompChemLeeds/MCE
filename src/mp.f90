MODULE mp

  use globvars

!*************************************************************************************************!
!*
!*         Morse Potential Module
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the Morse Potential Parameters
!*      2) Calculating the real and imaginary parts of the m dimensional initial wavefunction zinit
!*      3) Calculating a single element of the Hamiltonian matrix
!*      4) Calculating the derivative of the hamiltonian
!*      
!*************************************************************************************************!

contains

!--------------------------------------------------------------------------------------------------

  subroutine readparams_mp
    implicit none
    character(LEN=100)::LINE
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=128, file='inham.dat', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening inham.dat file'
      errorflag = 1
      return
    end if

    read(128,*,iostat=ierr)LINE

    do while (ierr==0)
      if(LINE=='MPw') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,freq_mp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading frequency value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='MPmass') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,mass_mp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading mass value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='MPDissEn') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,dissen_mp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Dissasociation Energy value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='MPWellParam') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,a0_mp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading well shaping parameter value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='MPupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='MPdownnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,lowlimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading lower limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      end if
      read (128,*,iostat=ierr)LINE
    end do

    close (128)

    if (n.ne.6) then
      write(0,"(a)") "Not all required variables read in readparams_mp subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_mp

!--------------------------------------------------------------------------------------------------

  subroutine genzinit_mp(mup, muq)   !   Level 1 Subroutine

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

    if (errorflag .ne. 0) return

    muq(1:ndim) = 5.0d0*sigq
    mup(1:ndim) = 0.0d0*sigp

    return

  end subroutine genzinit_mp

!--------------------------------------------------------------------------------------------------

  subroutine Hij_mp(H,z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    integer :: m
    complex(kind=8), dimension (:), allocatable :: Htemp
    complex(kind=8) :: zcomb

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Morse Potential"
      errorflag = 1
      return
    end if

    do m=1,ndim
      zcomb    = (dconjg(z1(m))+z2(m))/sqrt(2.0*gam)
      Htemp(m) = (-1.0d0/(4.0d0*mass_mp*mass_mp*freq_mp))*(dconjg(z1(m))**2+z2(m)**2-2*dconjg(z1(m))*z2(m)-1)
      Htemp(m) = Htemp(m) + dissen_mp*(1 + exp(a0_mp*(a0_mp-2.0d0*zcomb)))
      Htemp(m) = Htemp(m) - dissen_mp*2*exp((0.25d0*a0_mp**2.0d0)-(2.0d0*a0_mp*zcomb))
    end do

    H(1,1) = sum(Htemp(1:ndim))

    return   

  end subroutine Hij_mp

!--------------------------------------------------------------------------------------------------

  function dh_dz_mp(z)

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_mp
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8) :: dhdztmp, zcomb
    real(kind=8) :: fact
    integer :: m

    if (errorflag .ne. 0) return

    fact =2.0d0*dissen_mp*a0_mp/(sqrt(2.0d0*gam))

    do m=1,ndim   
      zcomb = (dconjg(z(m))+z(m))/sqrt(2.0*gam)  
      dhdztmp = (-1.0d0/(2.0d0*mass_mp*mass_mp*freq_mp))*(dconjg(z(m))-z(m))
      dhdztmp = dhdztmp - fact*exp(a0_mp*(a0_mp-2.0d0*zcomb))
      dhdztmp = dhdztmp + fact*exp((0.25d0*a0_mp**2.0d0)-(2.0d0*a0_mp*zcomb))
      dh_dz_mp (1,1,m) = dhdztmp
    end do

    return

  end function dh_dz_mp

!*************************************************************************************************!

end module mp
