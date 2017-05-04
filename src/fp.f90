MODULE fp

  use globvars

!***********************************************************************************!
!*
!*         Free Particle Module
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the Free Particle Parameters
!*      2) Calculating the real and imaginary parts of the m dimensional initial 
!*               wavefunction zinit
!*      3) Calculating a single element of the Hamiltonian matrix
!*      4) Calculating the derivative of the hamiltonian
!*      
!***********************************************************************************!

contains

!------------------------------------------------------------------------------------

  subroutine readparams_fp
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
      if(LINE=='FPmass') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,mass_fp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Mass value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='FPupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
          end if
          n = n+1
      else if(LINE=='FPdownnorm') then
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

    if (n.ne.3) then
      write(0,"(a)") "Not all required variables read in readparams_fp subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_fp

!------------------------------------------------------------------------------------

  subroutine genzinit_fp(mup, muq)   !   Level 1 Subroutine

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

    if (errorflag .ne. 0) return

    muq(1:ndim) = 5.0d0*sigp
    mup(1:ndim) = 0.0d0*sigq

    return

  end subroutine genzinit_fp

!------------------------------------------------------------------------------------

  subroutine Hij_fp(H,z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Harmonic Potential"
      errorflag = 1
      return
    end if

    H(1,1) = (-1.0d0/(4.0d0*gam*mass_fp*hbar))*&
           sum(dconjg(z1(1:ndim))**2+z2(1:ndim)**2-2*dconjg(z1(1:ndim))*z2(1:ndim)-1)

    return   

  end subroutine Hij_fp

!------------------------------------------------------------------------------------

  function dh_dz_fp(z)

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_fp
    complex(kind=8),dimension(:),intent(in)::z
    integer :: m

    if (errorflag .ne. 0) return

    do m=1,ndim    
      dh_dz_fp(1,1,m) = (-1.0d0/(2.0d0*gam*mass_fp*hbar))*(dconjg(z(m))-z(m))
    end do

    return

  end function dh_dz_fp

!***********************************************************************************!

end module fp
