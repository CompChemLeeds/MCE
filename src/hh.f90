MODULE hh

  use globvars

!***********************************************************************************!
!*
!*        Hennon-Heiles Potential Module
!*          
!*   Contains subroutines for:
!*
!*      1) Reading the Hennon-Heiles Potential Parameters
!*      2) Calculating the real and imaginary parts of the m dimensional initial 
!*          wavefunction zinit
!*      3) Calculating a single element of the Hamiltonian matrix
!*      4) Calculating the derivative of the hamiltonian
!*      
!***********************************************************************************!

contains

!------------------------------------------------------------------------------------

  subroutine readparams_hh
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
      if(LINE=='HHCoupling') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,lambda_hh
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading coupling constant value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='HHupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='HHdownnorm') then
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
      write(0,"(a)") "Not all required variables read in readparams_hh subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_hh

!------------------------------------------------------------------------------------

  subroutine genzinit_hh(mup, muq)   !   Level 1 Subroutine

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

    if (errorflag .ne. 0) return

    muq(1:ndim) = 2.0d0*sigq
    mup(1:ndim) = 0.0d0*sigp

    return

  end subroutine genzinit_hh

!------------------------------------------------------------------------------------

  subroutine Hij_hh(H,z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    complex(kind=8), dimension (:), allocatable :: z1c, Htemp
    real (kind=8) :: prefac
    integer :: m, ierr

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Harmonic Potential"
      errorflag = 1
      return
    end if

    allocate (z1c(ndim), stat=ierr)
    if (ierr==0) allocate (Htemp(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating z1c or Htemp in Hij_hh"
      errorflag=1
    end if

    prefac = lambda_hh/(2.0d0*dsqrt(2.0d0))

    do m=1,ndim
      z1c(m)=dconjg(z1(m))
    end do

    do m=1,ndim
      Htemp(m) = (z1c(m)*z2(m))+0.5d0
    end do

    do m=1,(ndim-1)
      Htemp(m) = Htemp(m) + prefac*((z1c(m)**2.0)+(z2(m)**2.0)&
                                          +(2.0*z1c(m)*z2(m))+1)*(z1c(m+1)+z2(m+1))
      Htemp(m) = Htemp(m) - prefac*((1.0/3.0)*((z1c(m+1)**3.0)+(z2(m+1)**3.0))&
                                                    +(3.0*(z1c(m+1)**2.0)*z2(m+1)))
      Htemp(m) = Htemp(m) - prefac*((3.0*z1c(m+1)*(z2(m+1)**2.0))+z1c(m+1)+z2(m+1))
    end do

    H(1,1) = sum(Htemp)

    deallocate (z1c, stat=ierr)
    if (ierr==0) deallocate (Htemp, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating z1c or Htemp in Hij_hh"
      errorflag=1
    end if

    return 

  end subroutine Hij_hh

!------------------------------------------------------------------------------------

  function dh_dz_hh(z)

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_hh
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8),dimension(:),allocatable::zc
    real(kind=8) :: fact
    integer :: m, ierr

    if (errorflag .ne. 0) return

    allocate (zc(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating zc array in dh_dz_hh function"
      errorflag = 1
      return
    end if

    fact = lambda_hh/(2.0*dsqrt(2.0d0))

    do m=1,ndim
      zc(m) = dconjg(z(m))
    end do

    m=ndim

    dh_dz_hh (1,1,1) = z(1) + (2.0*fact*(zc(1)+z(1))*(zc(2)+z(2)))
    dh_dz_hh (1,1,m) = z(m) + fact*((zc(m-1)**2.0)+(z(m-1)**2.0)&
                     +(2.0*zc(m-1)*z(m-1))-(zc(m)**2.0)-(z(m)**2.0)-(2.0*zc(m)*z(m)))

    if (ndim.gt.2) then
      do m=2,ndim-1 
        dh_dz_hh (1,1,m) = z(m) + fact*((zc(m-1)**2.0)+(z(m-1)**2.0)&
                         +(2.0*zc(m-1)*z(m-1))-(zc(m)**2.0)-(z(m)**2.0)&
                         -(2.0*zc(m)*z(m)) + (2.0*(zc(m)+z(m))*(zc(m+1)+z(m+1)))) 
      end do
    end if

    deallocate (zc, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating zc array in dh_dz_hh function"
      errorflag = 1
      return
    end if

    return

  end function dh_dz_hh

!***********************************************************************************!

end module hh
