MODULE hp

   use globvars

!***********************************************************************************!
!*
!*  Harmonic Potential Module
!*           
!*  Contains subroutines for:
!*
!*  1) Reading the Harmonic Potential Parameters
!*  2) Calculating the real and imaginary parts of the m dimensional initial 
!*      wavefunction zinit
!*  3) Calculating a single element of the Hamiltonian matrix
!*  4) Calculating the derivative of the hamiltonian
!*      
!***********************************************************************************!
  
contains

!------------------------------------------------------------------------------------

   subroutine readparams_hp
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
      if(LINE=='HPw') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,freq_hp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Frequency value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='HPupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='HPdownnorm') then
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
      write(0,"(a)") "Not all required variables read in readparams_hp subroutine"
      errorflag = 1
      return
    end if

    return

   end subroutine readparams_hp

!------------------------------------------------------------------------------------

   subroutine genzinit_hp(mup, muq) ! Level 1 Subroutine

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

    if (errorflag .ne. 0) return

    muq(1:ndim) = 10.0d0*sigp
    mup(1:ndim) = 0.0d0*sigq

    return

   end subroutine genzinit_hp

!------------------------------------------------------------------------------------

   subroutine Hij_hp(H,z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Harmonic Potential"
      errorflag = 1
      return
    end if

    H(1,1) = hbar*freq_hp*sum(dconjg(z1(1:ndim))*z2(1:ndim))+(0.5*ndim)

    return 

   end subroutine Hij_hp

!------------------------------------------------------------------------------------

   function dh_dz_hp(z)

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_hp
    complex(kind=8),dimension(:),intent(in)::z
    integer :: m

    if (errorflag .ne. 0) return

    do m=1,ndim
      dh_dz_hp(1,1,m) = hbar*freq_hp*z(m)
    end do

    return

   end function dh_dz_hp
   
!------------------------------------------------------------------------------------

  function disp_hp (bs)
  
    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8), dimension(:), allocatable :: D, Dc, d_small, dc_small, zk, zj, rho
    complex(kind=8) :: disp_hp, rho2, rhosum, ovrlp, term1, term2
    real(kind=8), dimension (:), allocatable :: s 
    integer :: j, k, m, ierr
    
    if (errorflag .ne. 0) return

    ierr = 0

    disp_hp = (0.0d0,0.0d0)

    if (ndim.ne.1) then
      write(0,"(a)") "Error! ndim should be 1 currently."
      errorflag=1
      return
    end if

    allocate (D(size(bs)), stat = ierr)
    if (ierr==0) allocate (Dc(size(bs)), stat=ierr)
    if (ierr==0) allocate (d_small(size(bs)), stat=ierr)
    if (ierr==0) allocate (dc_small(size(bs)), stat=ierr)
    if (ierr==0) allocate (zk(ndim), stat=ierr)
    if (ierr==0) allocate (zj(ndim), stat=ierr)
    if (ierr==0) allocate (rho(ndim), stat=ierr)
    if (ierr==0) allocate (s(size(bs)), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating basis set variable arrays for dispersion calculation"
      errorflag = 1
      return
    end if

    do k=1,size(bs)
      D(k)=bs(k)%D_big
      Dc(k)=dconjg(D(k))
      d_small(k)=bs(k)%d_pes(1)
      dc_small(k)=dconjg(d_small(k))
      s(k)=bs(k)%s_pes(1)
    end do         
    
    do k=1,size(bs)
      do j=1,size(bs)
        rho2 = (0.0d0,0.0d0)
        do m=1,ndim
          zk(m) = bs(k)%z(m)
          zj(m) = bs(j)%z(m)
          rho(m) = (dconjg(zj(m))+zk(m))/sqrt(2.0*gam)
        end do
        rho2 = sum(rho(1:ndim)**2.0d0)
        rhosum = sum(rho(1:ndim))
        ovrlp = product(cdexp((dconjg(zj(1:ndim))*zk(1:ndim))&
              -(0.5d0*dconjg(zj(1:ndim))*zj(1:ndim))&
              -(0.5d0*dconjg(zk(1:ndim))*zk(1:ndim))))
        term1 = Dc(j)*D(k)*dc_small(j)*d_small(k)*cdexp(i*(s(k)-s(j)))* ovrlp * rhosum
        term2 = Dc(j)*D(k)*dc_small(j)*d_small(k)*cdexp(i*(s(k)-s(j)))* ovrlp * (1./2.*gam + rho2)
        disp_hp = disp_hp + sqrt((term1 ** 2.) - term2)
      end do
    end do
    
    deallocate (D, stat = ierr)
    if (ierr==0) deallocate (Dc, stat=ierr)
    if (ierr==0) deallocate (d_small, stat=ierr)
    if (ierr==0) deallocate (dc_small, stat=ierr)
    if (ierr==0) deallocate (zk, stat=ierr)
    if (ierr==0) deallocate (zj, stat=ierr)
    if (ierr==0) deallocate (rho, stat=ierr)
    if (ierr==0) deallocate (s, stat=ierr)
    if (ierr/=0) then
      write(0,"(a,a)") "Error deallocating basis set variable arrays ",&
                 "for dipole acceleration calculation"   
      errorflag = 1   
      return
    end if
  
  end function disp_hp

!***********************************************************************************!

end module hp
