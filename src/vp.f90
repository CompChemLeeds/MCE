MODULE vp

  use globvars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                   !!!!!!!!!
!!!!!  NOTE: THIS MODULE HAS ONLY BEEN SET UP FOR 1D    !!!!!!!!!
!!!!!                                                   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!*************************************************************************************************!
!                                                                                                 !
!      This Module contains the subroutines and functions for the Spin Boson model.               !
!                                                                                                 !
!      It can be used as a template for other models, as all other subroutines and functions      !
!      in the program are very general, allowing different Hamiltonians or numbers of PESs        !
!      to be used.                                                                                !
!                                                                                                 !
!      The subroutines / functions which must be included in any future                           !
!      adaptations of this module are:                                                            !
!                                                                                                 !
!         1) ***subroutine readparams*** This reads the system specific parameters                !
!         2) ***subroutine genzinit*** This creates the initial CS                                !
!         3) ***function Hij(z1,z2)*** This builds a npes x npes matrix which represent a single  !
!                  basis function combination for the Hamiltonian                                 !
!         4) ***function dh_dz(z)***   This calculates the derivative of the Ehrenfest            !
!                  Hamiltonian for use in calculating the time derivative for z                   !
!                                                                                                 !
!       All of the subroutines and functions which are dependent upon the above should be         !
!       included in this module also. The subroutines should be referenced in the redirect        !
!       switchboard module and identified by a two character code, ie spin boson = sb             !
!                                                                                                 !
!*************************************************************************************************!

!*************************************************************************************************!
!*
!*         Spin Boson Module
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the Spin Boson Parameters
!*      2) Calculating the array for the frequency distribution wm
!*      3) Calculating the array for the amplitudes Cm 
!*      4) Calculating the array for the width distribution sig 
!*      5) Calculating the real and imaginary parts of the m dimensional initial wavefunction zinit
!*      6) Calculating a single 2 x 2 matrix which is equivalent to a single element of the 
!*         Hamiltonian matrix
!*      7) Calculating the Bath Hamiltonian for a pair of basis functions
!*      8) Calculating the Coupling Hamiltonian for a pair of basis functions
!*      9) Combining the elements of the differential of the hamiltonan for each combination of 2
!*         PESs with respect to z
!*     10) Calculating the differential of the hamiltonian for the PES combination (1,1)
!*     11) Calculating the differential of the hamiltonian for the PES combination (1,2)
!*     12) Calculating the differential of the hamiltonian for the PES combination (2,1)
!*     13) Calculating the differential of the hamiltonian for the PES combination (2,2)
!*      
!*************************************************************************************************!

contains

!--------------------------------------------------------------------------------------------------

  subroutine readparams_vp
  
    !   This subroutine reads the system specific parameters from the inham.dat file
  
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

      if(LINE=='VP_y_00') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,y_00_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading S_0 offset value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='VP_A_0') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,A_0_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading S_0 amplitude value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='VP_b_0') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,b_0_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading S_0 parabolic parameter value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='VPwidth') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,width_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading S_0 Gaussian width value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='VP_A_NAC') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,A_NAC_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Non-adiabatic Coupling amplitude value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='VP_y_10') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,y_10_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading S_1 offset value"
          errorflag = 1
          return
        end if
        n = n+1          
      else if (LINE=='VP_b_1') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,b_1_vp
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading S_1 parabolic parameter value"
          errorflag = 1
          return
        end if
        n = n+1                
      else if(LINE=='VPupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='VPdownnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,lowlimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading lower limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      end if

      read(128,*,iostat=ierr) LINE

    end do

    close (128)

    if (n.ne.9) then
      write(0,"(a)") "Not all required variables read in readparams_vp subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_vp

!--------------------------------------------------------------------------------------------------

  subroutine genzinit_vp(mup, muq)   !   Level 1 Subroutine

    ! Generates the initial Gaussian, used as a central point for the larger basis set (which is calculated in a different function)
    
    implicit none

    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq
    integer :: m

    if (errorflag .ne. 0) return
    
    if ((.not.allocated(mup)).or.(.not.allocated(muq)).or.(size(mup).ne.size(muq)) &
           .or.(size(mup).ne.ndim)) then
      write(0,"(a)") 'Allocation error in mup or muq'
      errorflag = 1
      return
    end if

    do m=1,ndim
      muq(m) = 5.0d0
      mup(m) = 0.0d0
    end do

    return

  end subroutine genzinit_vp

!--------------------------------------------------------------------------------------------------

  subroutine Hij_vp(H,z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    complex(kind=8):: HS0, HS1, HNAC
    real(kind=8) :: chk
    integer :: ierr

    if (errorflag .ne. 0) return

    ierr = 0

    HS0 = HS0_vp(z1,z2)
    HS1 = HS1_vp(z1,z2)
    HNAC = HNAC_vp(z1,z2)

    chk = sum((abs(z1(1:ndim)-z2(1:ndim)))**2)

    if (chk.le.20000) then
      H(1,1) = HS0
      H(1,2) = HNAC
      H(2,1) = HNAC
      H(2,2) = HS1
    else
      H = (0.0d0, 0.0d0)
    end if

    return   

  end subroutine Hij_vp

!--------------------------------------------------------------------------------------------------

  function HS0_vp(z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8)::HS0_vp, rho2, rhosum, Htemp
    complex (kind=8), dimension (:), allocatable :: rho, z1c
    real(kind=8) :: fact
    integer :: m

    if (errorflag .ne. 0) return

    HS0_vp = (0.0d0, 0.0d0)
    Htemp = (0.0d0,0.0d0)
    
    allocate (rho(size(z1)))
    allocate (z1c(size(z1)))

    do m=1,ndim
      z1c(m) = dconjg(z1(m))
      rho(m)=((z1c(m)+z2(m))/2*gam)
    end do
    
    rho2 = sum(rho(1:ndim)**2.0d0)
    rhosum = sum(rho(1:ndim))

    do m=1,ndim
      Htemp = Htemp - ((hbar**2)*gam/(4.0d0))*(z1c(m)**2.0d0+z2(m)**2.0d0-2.0d0*z1c(m)*z2(m)-1.0d0)&
                            + (y_00_vp - A_0_vp) + (b_0_vp/(2.*gam)+b_0_vp*rho(m)**2.)
    end do
       
    fact = (A_0_vp*width_vp/sqrt(width_vp**2.+(0.5/gam)))**dble(ndim)
    
!    Htemp = sqrt(2.0d0*b_0_vp)*sum(z1c(1:ndim)*z2(1:ndim))+(0.5*ndim)
    
    Htemp = Htemp+product(fact*cdexp(-gam*rho(1:ndim)**2./(2.*gam*width_vp+1.)))
    
    HS0_vp = Htemp   
    
    deallocate(rho,z1c) 

    return

  end function HS0_vp
    
!--------------------------------------------------------------------------------------------------

  function HS1_vp(z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8)::HS1_vp
    complex (kind=8), dimension (:), allocatable :: rho, z1c
    integer :: m

    if (errorflag .ne. 0) return

    allocate (rho(size(z1)))
    allocate (z1c(size(z1)))

    HS1_vp = (0.0d0, 0.0d0)

    do m=1,ndim
      z1c(m)=dconjg(z1(m))
      rho(m)=((z1c(m)+z2(m))/2*gam)
    end do
       
    do m=1,ndim
      HS1_vp = HS1_vp - ((hbar**2)*gam/(4.0d0))*(z1c(m)**2.0d0+z2(m)**2.0d0-2.0d0*z1c(m)*z2(m)-1.0d0)&
                            + (y_10_vp) + (b_1_vp/(2.*gam)+b_1_vp*rho(m)**2.)
    end do
   
    deallocate(rho,z1c)
   
!    HS1_vp = sqrt(2.0d0*b_1_vp)*(sum(dconjg(z1(1:ndim))*z2(1:ndim))+(0.5*ndim))+(y_00_vp*ndim)

    return

  end function HS1_vp

!--------------------------------------------------------------------------------------------------

  function HNAC_vp(z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8)::HNAC_vp, Htemp
    complex (kind=8), dimension (:), allocatable :: rho
    integer :: m

    if (errorflag .ne. 0) return
    
    allocate (rho(size(z1)))

    HNAC_vp = (0.0d0, 0.0d0)
    Htemp = (0.0d0,0.0d0)

    do m=1,ndim
      rho(m)=((dconjg(z1(m))+z2(m))/2*gam)
    end do
    
    Htemp = Htemp + product(sqrt((A_NAC_vp*gam/(gam+1.)))*cdexp(-dble(ndim)*gam*rho(1:ndim)**2./(gam+1.)))
    
    HNAC_vp = Htemp   
    
    deallocate (rho) 

    return

  end function HNAC_vp
  
!---------------------------------------------------------------------------------------------------  

  function dh_dz_vp(z)

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_vp
    complex(kind=8),dimension(:),intent(in)::z
    integer :: ierr

    if (errorflag .ne. 0) return

    if (npes.ne.2) then
      write(0,"(a)") "No. PES is not equal to 2, but Pseudo-Tortion 2 PES derivative called"
      errorflag = 1
      return
    end if

    dh_dz_vp(1,1,:) = dhdz_vp_11(z)
    dh_dz_vp(1,2,:) = dhdz_vp_12(z)
    dh_dz_vp(2,1,:) = dhdz_vp_21(z)
    dh_dz_vp(2,2,:) = dhdz_vp_22(z)

    return

  end function dh_dz_vp

!--------------------------------------------------------------------------------------------------

  function dhdz_vp_11(z)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_vp_11
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8), dimension (:), allocatable :: zc, rho
    real(kind=8) :: fact
    integer::m

    if (errorflag .ne. 0) return

    allocate (zc(size(z)))
    allocate (rho(size(z)))
    zc(1:ndim)=dconjg(z(1:ndim))

    do m=1,ndim
      rho(m) = (zc(m)+z(m))/sqrt(2.0*gam)
    end do
    
    fact = (2.*A_0_vp*width_vp/sqrt(4.*width_vp**2.+(2./gam)))/(2.*gam*width_vp**2.+1.)
    
    do m=1,ndim
      dhdz_vp_11(m) = (0.0d0,0.0d0)
      dhdz_vp_11(m) = dhdz_vp_11(m) - (gam/2.0d0)*(zc(m)-z(m)) 
      dhdz_vp_11(m) = dhdz_vp_11(m) + (b_0_vp*rho(m)/gam)-(fact*rho(m)*cdexp(-1.*gam*rho(m)**2./(2.*gam*width_vp**2.+1.)))
!      dhdz_vp_11(m) = dhdz_vp_11(m) + sqrt(2.0d0*b_0_vp)*z(m)-(fact*rho(m)*cdexp(-1.*gam*rho(m)**2./(2.*gam*width_vp**2.+1.)))
    end do

    deallocate(zc)
    deallocate(rho)

    return

  end function dhdz_vp_11

!--------------------------------------------------------------------------------------------------

  function dhdz_vp_12(z)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_vp_12
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8), dimension (:), allocatable :: zc, rho
    integer::m

    if (errorflag .ne. 0) return
    
    allocate (zc(size(z)))
    allocate (rho(size(z)))
    zc(1:ndim)=dconjg(z(1:ndim))

    do m=1,ndim
      rho(m) = (zc(m)+z(m))/sqrt(2.0*gam)
    end do

    do m=1,ndim
      dhdz_vp_12(m) = (0.0d0,0.0d0)
      dhdz_vp_12(m) = dhdz_vp_12(m)-(sqrt(A_NAC_vp**2.*gam/(gam+1.)))*rho(m)/(gam+1.)*cdexp(-1.*gam*rho(m)**2./(gam+1.))
    end do

    deallocate(zc)
    deallocate(rho)
    
    return

  end function dhdz_vp_12

!--------------------------------------------------------------------------------------------------

  function dhdz_vp_21(z)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_vp_21
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8), dimension (:), allocatable :: zc, rho
    integer::m

    if (errorflag .ne. 0) return
    
    allocate (zc(size(z)))
    allocate (rho(size(z)))
    zc(1:ndim)=dconjg(z(1:ndim))

    do m=1,ndim
      rho(m) = (zc(m)+z(m))/sqrt(2.0*gam)
    end do

    do m=1,ndim
      dhdz_vp_21(m) = (0.0d0,0.0d0)
      dhdz_vp_21(m) = dhdz_vp_21(m)-(sqrt(A_NAC_vp**2.*gam/(gam+1.)))*rho(m)/(gam+1.)*cdexp(-1.*gam*rho(m)**2./(gam+1.))
    end do

    deallocate(zc)
    deallocate(rho)
    
    return

  end function dhdz_vp_21

!--------------------------------------------------------------------------------------------------

  function dhdz_vp_22(z)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_vp_22
    complex(kind=8),dimension(:),intent(in)::z
    complex(kind=8), dimension (:), allocatable :: zc, rho
    integer::m

    if (errorflag .ne. 0) return
    
    allocate (zc(size(z)))
    allocate (rho(size(z)))
    zc(1:ndim)=dconjg(z(1:ndim))

    do m=1,ndim
      rho(m) = (zc(m)+z(m))/sqrt(2.0*gam)
    end do
    
    do m=1,ndim
      dhdz_vp_22(m) = (0.0d0,0.0d0)
      dhdz_vp_22(m) = dhdz_vp_22(m) - (gam/2.0d0)*(zc(m)-z(m)) + (b_1_vp*rho(m)/gam)     
!      dhdz_vp_22(m) = dhdz_vp_22(m) + sqrt(2.0d0*b_1_vp)*z(m)
    end do

    deallocate(zc)
    deallocate(rho)
    
    return

  end function dhdz_vp_22

!*************************************************************************************************!

  function disp_vp (bs)
  
    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8), dimension(:), allocatable :: D, Dc, d_small, dc_small, zk, zj, rho
    complex(kind=8) :: disp_vp, rho2, rhosum, ovrlp, term1, term2
    real(kind=8), dimension (:), allocatable :: s 
    integer :: j, k, m, ierr
    
    if (errorflag .ne. 0) return

    ierr = 0

    disp_vp = (0.0d0,0.0d0)

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
        disp_vp = disp_vp + sqrt((term1 ** 2.) - term2)
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
  
  end function disp_vp

!*************************************************************************************************!

  function gauss_random_vp (width, muq, mup)
  
    implicit none
    complex(kind=8) :: gauss_random_vp
    real(kind=8), intent(in) :: width, mup, muq
    integer :: iset
    real(kind=8) :: fac,rsq,v1,v2,ran_1,xz,yz

    iset = 0
    
    do while (iset .lt. 2)      
      CALL RANDOM_NUMBER(ran_1) 
      v1=2.0d0*ran_1-1.0d0
      CALL RANDOM_NUMBER(ran_1)     
      v2=2.0d0*ran_1-1.0d0
      rsq=v1**2.0d0+v2**2.0d0
      if (rsq.ge.1.0d0.or.rsq.eq.0.0d0) cycle
      fac=sqrt(-2.0d0*log(rsq)/rsq)
      if (iset == 0) then
        xz=v2*fac*sigq
      else if (iset == 1) then
        yz=v2*fac*sigp
      end if    
      iset=iset + 1
    end do
    
    xz = muq + (xz * width)
    yz = mup + (yz * width)
    
    gauss_random_vp = cmplx(xz,yz,kind=8)
    
    return
    
  end function gauss_random_vp   

!*************************************************************************************************!

end module vp
