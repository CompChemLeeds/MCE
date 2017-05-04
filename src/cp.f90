MODULE cp

  use double_complex_erf
  use globvars

!***********************************************************************************!
!*
!*         Coulomb Potential Module (For use with a strong laser field to show HHG)
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the Coulomb Potential Parameters
!*      2) Calculating the real and imaginary parts of the m dimensional initial 
!*          wavefunction zinit
!*      3) Calculating a single element of the Hamiltonian matrix
!*      4) Calculating the derivative of the hamiltonian
!*      
!***********************************************************************************!

contains

!------------------------------------------------------------------------------------

  subroutine readparams_cp

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
      if(LINE=='CPmass') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,mass_cp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading mass value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='CPfrequency') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,freq_cp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading frequency value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='CPIntensity') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,inten_cp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading laser intensity value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='CPRc') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,Rc_cp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading distance-to-centre value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='CPupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='CPdownnorm') then
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

    if (n.ne.5) then
      write(0,"(a)") "Not all required variables read in readparams_cp subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_cp

!------------------------------------------------------------------------------------

  subroutine genzinit_cp(mup, muq)

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq
    
    integer :: m

    if (errorflag .ne. 0) return
    
    do m=1,ndim
!      if (m==1) then
!        muq(m) = inten_cp/(freq_cp**2)*sigq
!      else
        muq(m) = 0.0d0*sigq
!      end if
      mup(m) = 0.0d0*sigp
    end do    

    return

  end subroutine genzinit_cp

!------------------------------------------------------------------------------------

  subroutine Hij_cp(H,z1,z2,t)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    real(kind=8), intent (in) :: t
    integer :: m, ierr
    complex(kind=8), dimension (:), allocatable :: Htemp, rho
    complex(kind=8) :: rho2, z1c 

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Inverse Gaussian"
      errorflag = 1
      return
    end if

    if (ndim.ne.3) then
      write(0,"(a)") "Error! Coulomb potential is only valid in 3D"
      errorflag = 1
      return
    end if

    allocate (Htemp(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating Htemp in Hij_cp"
      errorflag=1
    end if

    allocate(rho(size(z1)))
    do m=1,ndim
      rho(m)=((dconjg(z1(m))+z2(m))/sqrt(2.0d0*gam))-Rc_cp
    end do
    rho2 = sum(rho(1:ndim)**2.0d0)

    do m=1,ndim
      z1c = dconjg(z1(m))
      Htemp(m) = (0.0d0,0.0d0)
      Htemp(m) = Htemp(m) - (1.0d0/4.0d0)*&
                       (z1c**2.0d0+z2(m)**2.0d0-2.0d0*z1c*z2(m)-1.0d0) !free particle
      if (abs(rho2) < tiny(1.0d0)) then             !checks that not near singularity
        Htemp(m) = Htemp(m) + 2.0d0*sqrt(gam)/(sqrtpi)   !empirically determined limit
      else
        Htemp(m) = Htemp(m) + (c_error_func(sqrt(gam*rho2)/sqrt(rho2))) !Atom. Pot.
      end if
      if (m==1)Htemp(m) = Htemp(m) + inten_cp*rho(m)*dcos(freq_cp*t)    !laser field
    end do

    H(1,1) = sum(Htemp(1:ndim))

    deallocate(Htemp, rho)

    return   

  end subroutine Hij_cp

!------------------------------------------------------------------------------------

  function dh_dz_cp(z,t) 

! This is laid out so as to allow easy conversion to calculating off diagonal elements

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_cp
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8), dimension(:), allocatable :: rho
    real(kind=8), intent (in) :: t
    complex(kind=8) :: dhdztmp, rho2, arho, datpotdz
    real(kind=8) :: rt2g
    integer :: m

    if (errorflag .ne. 0) return

    rt2g = sqrt(2.0d0*gam)
    allocate(rho(size(z)))
    do m=1,ndim
      rho(m)=(dconjg(z(m))+z(m))/sqrt(2.0d0*gam)-Rc_cp
    end do
    rho2 = sum(rho(1:ndim)**2.0d0)
    arho = sqrt(rho2)

    datpotdz = (2.0d0*sqrt(gam)*cdexp(-gam*rho2))/(sqrtpi*rho2)-&
                                        (c_error_func(sqrt(gam)*arho)/(arho**3.0))

    do m=1,ndim    
      dhdztmp = (0.0d0,0.0d0) 
      dhdztmp = dhdztmp - (1.0/2.0)*(dconjg(z(m))-z(m))            ! free particle
      if (abs(rho2) .ge. tiny(1.0d0)) then             !checks that not near singularity
        if (datpotdz /= datpotdz) then  !check for NaN or Inf from the exponential
          continue                     !take no action
        else
          dhdztmp = dhdztmp + (1.0d0/(rt2g))*rho(m)*datpotdz  !atomic potential
        end if
      end if
      if (m==1) dhdztmp = dhdztmp + (inten_cp/rt2g)*dcos(freq_cp*t)  ! laser field
      dh_dz_cp (1,1,m) = dhdztmp
    end do

    return

  end function dh_dz_cp

!------------------------------------------------------------------------------------

    function dipole_cp(bs)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8) :: dipole_cp, ovrlp, rho2, arho, datpotdz
    complex(kind=8), dimension(:), allocatable :: D, Dc, zk, zj, rho
    real(kind=8), dimension (:), allocatable :: s
    integer::k,j,m,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    dipole_cp = (0.0d0,0.0d0)

    if ((ndim.ne.1).and.(ndim.ne.3)) then
      write(0,"(a)") "Error! ndim should be 1 or 3 but is not."
      errorflag=1
      return
    end if

    if ((npes.ne.1).or.(in_pes.ne.1)) then
      write(0,"(a)") "Error! npes and in_pes should be 1 but are not."
      errorflag=1
      return
    end if    

    allocate (D(size(bs)), stat = ierr)
    if (ierr==0) allocate (Dc(size(bs)), stat=ierr)
    if (ierr==0) allocate (zk(ndim), stat=ierr)
    if (ierr==0) allocate (zj(ndim), stat=ierr)
    if (ierr==0) allocate (rho(ndim), stat=ierr)
    if (ierr==0) allocate (s(size(bs)), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating basis set variable arrays for dipole calculation"
      errorflag = 1
      return
    end if

    do k=1,size(bs)
      D(k)=bs(k)%D_big
      Dc(k)=dconjg(D(k))
      s(k)=bs(k)%s_pes(1)
    end do     

    do k=1,size(bs)
      do j=1,size(bs)
        do m=1,ndim
          zk(m) = bs(k)%z(m)
          zj(m) = bs(j)%z(m)
          rho(m) = (dconjg(zj(m))+zk(m))/sqrt(2.0*gam)-Rc_cp
        end do
        rho2 = sum(rho(1:ndim)**2.0d0)
        arho = sqrt(rho2)
        datpotdz = (2.0d0*sqrt(gam)*cdexp(-gam*rho2))/(sqrtpi*rho2)-&
                                           (c_error_func(sqrt(gam)*arho)/(arho**3.0))
        ovrlp = product(cdexp((dconjg(zj(1:ndim))*zk(1:ndim))&
                              -(0.5d0*dconjg(zj(1:ndim))*zj(1:ndim))&
                              -(0.5d0*dconjg(zk(1:ndim))*zk(1:ndim))))
        if (abs(arho) > tiny(1.0d0)) then           !checks that not near singularity
          if (datpotdz /= datpotdz) then  !check for NaN or Inf from the exponential
            continue                     !take no action
          else
            do m=1,ndim
              dipole_cp = dipole_cp+(-1.0d0*Dc(j)*D(k)*cdexp(i*(s(k)-s(j)))* ovrlp*&
                 (1.0d0/(sqrt(gam)*sqrtpi))*rho(m)*datpotdz)
            end do
          end if
        end if
      end do
    end do

    deallocate (D, stat = ierr)
    if (ierr==0) deallocate (Dc, stat=ierr)
    if (ierr==0) deallocate (zk, stat=ierr)
    if (ierr==0) deallocate (zj, stat=ierr)
    if (ierr==0) deallocate (rho, stat=ierr)
    if (ierr==0) deallocate (s, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating basis set variable arrays for dipole calculation"   
      errorflag = 1   
      return
    end if

    return

  end function dipole_cp

!------------------------------------------------------------------------------------

  function c_error_func(z) ! wrapper for dcerf to filter out NaNs. From A. Kirrander

    implicit none
    complex(kind=8)    :: c_error_func ! output
    complex(kind=8)    :: z            ! input
    double precision   :: out(2)
    integer, parameter :: flag=0       ! 0=erfz, 1=complementary erfz

    call dcerf (flag, (/dble(z), aimag(z)/), out) 
    c_error_func = out(1) + i*out(2)

    if (c_error_func /= c_error_func ) then ! check for overflow (inf,nan) from dcerf
      c_error_func = (0.d0,0.d0)
    end if

  end function c_error_func

!************************************************************************************!

end module cp

