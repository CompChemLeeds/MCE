MODULE dl

  use globvars

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
!*         Spin Boson Module with drude lorentz spectral density
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

  subroutine readparams_dl
  
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

      if(LINE=='DLDelta') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,delta_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Delta value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='DLEps') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,eps_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Epsilon value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='DLw') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,wc_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wc value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='DLlambda') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,lambda_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading lambda parameter value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='DLwmax') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,wmax_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading wmax value"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='DLBeta') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,beta_dl
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading beta value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='DLupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='DLdownnorm') then
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

!    if (wmax_dl .ne. 5.0d0 * wc_dl) then
!      wmax_dl = wc_dl * 5.0d0
!    end if

    close (128)

    if (n.ne.8) then
      write(0,"(a)") "Not all required variables read in readparams_dl subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_dl

!--------------------------------------------------------------------------------------------------

  subroutine genwm_dl(wm_dl)   !   Level 2 Subroutine

    !   Generates the array of discretised frequencies over M degrees of freedom

    implicit none
    integer::m, ierr
    real(kind=8), dimension(:), allocatable, intent(inout) :: wm_dl

    if (errorflag .ne. 0) return
    ierr = 0

    if (ndim.le.0) then
      write(0,"(a)") "Number of dimensions not correctly read or stored"
      errorflag = 1
      return
    else if (.not.allocated(wm_dl)) then
      write(0,"(a)") "genwm subroutine called but wm not allocated"
      errorflag=1
      return
    else
      do m=1,size(wm_dl)
        wm_dl(m)=wc_dl*dtan((dble(m)/dble(ndim))*datan(wmax_dl/wc_dl))                           ! Drude-Lorentz (or Debye) spectral density
      end do
    end if

    return

  end subroutine genwm_dl

!--------------------------------------------------------------------------------------------------

  subroutine genCm_dl(Cm_dl, wm_dl)   !   Level 2 Subroutine
  
    ! Generates the array of coupling strengths, which are dependant upon the frequencies, also over M degrees of freedom

    implicit none
    real(kind=8) :: pi_rl
    integer::m, ierr
    real(kind=8), dimension(:), intent(in) :: wm_dl 
    real(kind=8), dimension(:), allocatable, intent(inout) :: Cm_dl  

    if (errorflag .ne. 0) return
    ierr = 0
    pi_rl = sqrtpi**2.0d0

    if (ndim.le.0) then
      write(0,"(a)") 'Number of dimensions not correctly read or stored'
      errorflag = 1
      return
    else if (.not.allocated(Cm_dl)) then
      write(0,"(a)") "genCm subroutine called but Cm not allocated"
      errorflag=1
      return
    else
      do m=1,size(Cm_dl)
        Cm_dl(m)=wm_dl(m)*sqrt((4.0d0*lambda_dl/(pi_rl*wc_dl))*datan(wmax_dl/wc_dl)/dble(ndim))
        ! ^^^ Drude-Lorentz (or Debye) spectral density.
      end do
    end if

    return

  end subroutine genCm_dl

!--------------------------------------------------------------------------------------------------

  subroutine gensig_dl(sig_dl, wm_dl)   !   Level 2 Subroutine

    ! Generates the distribution of states over the M degrees of freedom in the initial gaussian.

    implicit none
    real(kind=8), dimension (:), allocatable, intent (inout) :: sig_dl
    real(kind=8), dimension (:), intent (in) :: wm_dl
    integer::m, ierr

    if (errorflag .ne. 0) return
    ierr = 0

    if (ndim.le.0) then
      write(0,"(a)") 'Number of dimensions not correctly read or stored'
      errorflag = 1
      return
    else if (.not.allocated(sig_dl)) then
      write(0,"(a)") 'Allocation error in sig'
      errorflag = 1
      return
    else
      do m=1,size(sig_dl)
        sig_dl(m) = 1.0d0/sqrt(dexp(beta_dl*wm_dl(m))-1.0d0)
      end do
    end if

    return

  end subroutine gensig_dl

!--------------------------------------------------------------------------------------------------

  subroutine genzinit_dl(mup, muq)   !   Level 1 Subroutine

    ! Generates the initial Gaussian, used as a central point for the larger basis set (which is calculated in a different function)
    
    ! z_in is defined as being subject to the distribution rho(zin) = sigma_dl(m)*exp(-1.0*sigma_dl(m)*abs(zin(m)))

    implicit none
    integer::m, n, recalc, ierr
    real(kind=8) :: Ezin
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq
    real(kind=8), dimension(:), allocatable :: sig_dl, wm_dl, Cm_dl
    complex(kind=8), dimension (:), allocatable:: zin
    complex(kind=8),dimension(:,:), allocatable::H

    if (errorflag .ne. 0) return

    n = 0
    recalc = 1 ! As there is always the possibility of outliers when using a rng based system, the 
    ierr = 0   !ability to reclculate an unsuitable z_in is included

    allocate(zin(ndim), H(npes,npes), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of zin and H in genzinit"
      errorflag=1
      return
    end if  

    if (ndim.le.0) then
      write(0,"(a)") 'Number of dimensions not correctly read or stored'
      errorflag = 1
      return
    end if

    if ((.not.allocated(mup)).or.(.not.allocated(muq)).or.(size(mup).ne.size(muq)) &
           .or.(size(mup).ne.ndim)) then
      write(0,"(a)") 'Allocation error in mup or muq'
      errorflag = 1
      return
    end if

    allocate (sig_dl(ndim), stat = ierr)
    if (ierr==0) allocate (wm_dl(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_dl(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in sig or wm allocation in genzinit"
      errorflag=1
      return
    end if

    call genwm_dl(wm_dl)
    call gensig_dl(sig_dl, wm_dl)    

    do while (recalc == 1)
      do m=1,ndim
!        mup(m)=ZBQLNOR(mu,sig_dl(m)*sigp)
!        muq(m)=ZBQLNOR(mu,sig_dl(m)*sigq)
!        zin(m)=cmplx(muq(m),mup(m),kind=8)
        zin(m)=gauss_random_dl(sig_dl(m),mu,mu)
        muq(m)=dble(zin(m))
        mup(m)=dimag(zin(m))
      end do
      if (ECheck.eq."YES") then
        call Hij_dl(H,zin,zin)
        Ezin = dble(H(in_pes,in_pes))
        if ((Ezin.gt.Ebfmax).or.(Ezin.lt.Ebfmin)) then
          if (n.lt.Ntries) then
            n = n+1
            write(6,"(a)")"Initial Basis did not meet energy requirements. Recalculating..."
            write(6,"(3(a,e12.5))")"Ebfmax = ", Ebfmax, "Ebfmin = ", Ebfmin, "Ezin = ", Ezin
            recalc = 1
          else
            write(6,"(a)")"Initial basis recalculated ", n, "times but still outside acceptable region."
            write(0,"(a)")"Multiple recalculations due to energy checks resulted in an unsuitable basis"
            write(0,"(a)")"Terminating calculation"
            errorflag = 1
            recalc = 0
          end if
        else
          recalc = 0           
        end if
      else
        recalc = 0
      end if
    end do

    deallocate(zin, H, sig_dl, wm_dl, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of zin, H, wm or sig in genzinit"
      errorflag=1
      return
    end if 

    return

  end subroutine genzinit_dl

!--------------------------------------------------------------------------------------------------

  subroutine Hij_dl(H,z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    complex(kind=8):: Hbath, Hcoup
    real(kind=8), dimension(:), allocatable :: wm_dl, Cm_dl
    real(kind=8) :: chk
    integer :: ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate (wm_dl(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_dl(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in wm or Cm allocation in Hij"
      errorflag=1
      return
    end if

    call genwm_dl(wm_dl)
    call genCm_dl(Cm_dl, wm_dl)

    Hbath = Hb_dl(z1,z2,wm_dl)
    Hcoup = Hc_dl(z1,z2,wm_dl, Cm_dl)

    chk = sum((abs(z1(1:ndim)-z2(1:ndim)))**2)

    if (chk.le.20000) then
      H(1,1) = Hbath+eps_dl+Hcoup
      H(1,2) = cmplx(delta_dl,0.0d0,kind=8)!+Hcoup
      H(2,1) = cmplx(delta_dl,0.0d0,kind=8)!+Hcoup
      H(2,2) = Hbath-eps_dl-Hcoup
    else
      H = (0.0d0, 0.0d0)
    end if

    deallocate(wm_dl, Cm_dl, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of wm or Cm in Hij"
      errorflag=1
      return
    end if 

    return   

  end subroutine Hij_dl

!--------------------------------------------------------------------------------------------------

  function Hb_dl(z1,z2,wm_dl)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    real(kind=8), dimension(:), intent(in) :: wm_dl
    complex(kind=8)::Hb_dl
    integer :: m

    if (errorflag .ne. 0) return

    Hb_dl = (0.0d0, 0.0d0)

    do m=1,ndim
      Hb_dl = Hb_dl + (wm_dl(m)*(dconjg(z1(m))*z2(m)+0.5d0))
    end do

    return

  end function Hb_dl
    
!--------------------------------------------------------------------------------------------------

  function Hc_dl(z1,z2,wm_dl, Cm_dl)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    real(kind=8), dimension(:), intent(in) :: wm_dl, Cm_dl
    complex(kind=8)::Hc_dl
    integer :: m

    if (errorflag .ne. 0) return

    Hc_dl = (0.0d0, 0.0d0)

    do m=1,ndim
      Hc_dl = Hc_dl + ((Cm_dl(m)/sqrt(2.0d0*wm_dl(m)))*(dconjg(z1(m))+z2(m)))
    end do

    return

  end function Hc_dl

!--------------------------------------------------------------------------------------------------

  function dh_dz_dl(z)

    implicit none
    complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_dl
    complex(kind=8),dimension(:),intent(in)::z
    real(kind=8), dimension(:), allocatable :: wm_dl, Cm_dl
    integer :: ierr

    if (errorflag .ne. 0) return

    if (npes.ne.2) then
      write(0,"(a)") "No. PES is not equal to 2, but Spin Boson 2 PES derivative called"
      errorflag = 1
      return
    end if

    allocate (wm_dl(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_dl(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in wm or Cm allocation in dh_dz"
      errorflag=1
      return
    end if 

    call genwm_dl(wm_dl)
    call genCm_dl(Cm_dl, wm_dl)

    dh_dz_dl(1,1,:) = dhdz_dl_11(z,wm_dl,Cm_dl)
    dh_dz_dl(1,2,:) = dhdz_dl_12(wm_dl,Cm_dl)
    dh_dz_dl(2,1,:) = dhdz_dl_21(wm_dl,Cm_dl)
    dh_dz_dl(2,2,:) = dhdz_dl_22(z,wm_dl,Cm_dl)

!    dh_dz_dl(1,1,:) = dhdz_dl_11(z,wm_dl)
!    dh_dz_dl(1,2,:) = dhdz_dl_12(wm_dl,Cm_dl)
!    dh_dz_dl(2,1,:) = dhdz_dl_21(wm_dl,Cm_dl)
!    dh_dz_dl(2,2,:) = dhdz_dl_22(z,wm_dl)

    deallocate (wm_dl, Cm_dl, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in wm or Cm allocation in dh_dz"
      errorflag=1
      return
    end if 

    return

  end function dh_dz_dl

!--------------------------------------------------------------------------------------------------

  function dhdz_dl_11(z,wm_dl, Cm_dl)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_dl_11
    complex(kind=8),dimension(:),intent(in)::z
    real(kind=8),dimension(:),intent(in)::wm_dl, Cm_dl
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
!      dhdz_dl_11(m) = wm_dl(m)*z(m)
      dhdz_dl_11(m) = wm_dl(m)*z(m) + (Cm_dl(m)/sqrt(2.0d0*wm_dl(m)))
    end do

    return

  end function dhdz_dl_11

!--------------------------------------------------------------------------------------------------

  function dhdz_dl_12(wm_dl,Cm_dl)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_dl_12
    real(kind=8),dimension(:),intent(in)::wm_dl, Cm_dl
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
!      dhdz_dl_12 = (Cm_dl(m)/sqrt(2.0d0*wm_dl(m)))
      dhdz_dl_12 = (0.0d0,0.0d0)
    end do

    return

  end function dhdz_dl_12

!--------------------------------------------------------------------------------------------------

  function dhdz_dl_21(wm_dl,Cm_dl)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_dl_21
    real(kind=8),dimension(:),intent(in)::wm_dl, Cm_dl
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
!      dhdz_dl_21(m) = (Cm_dl(m)/sqrt(2.0d0*wm_dl(m)))
      dhdz_dl_21(m) = (0.0d0,0.0d0)
    end do

    return

  end function dhdz_dl_21

!--------------------------------------------------------------------------------------------------

  function dhdz_dl_22(z,wm_dl,Cm_dl)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_dl_22
    complex(kind=8),dimension(:),intent(in)::z
    real(kind=8),dimension(:),intent(in)::wm_dl, Cm_dl
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
      dhdz_dl_22(m) = wm_dl(m)*z(m) - (Cm_dl(m)/sqrt(2.0d0*wm_dl(m)))
!      dhdz_dl_22(m) = wm_dl(m)*z(m)
    end do

    return

  end function dhdz_dl_22

!*************************************************************************************************!

  function gauss_random_dl (width, muq, mup)
  
    implicit none
    complex(kind=8) :: gauss_random_dl
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
        xz=v2*fac/sqrt(2.0d0)
      else if (iset == 1) then
        yz=v2*fac/sqrt(2.0d0)
      end if    
      iset=iset + 1
    end do
    
    xz = muq + (xz * width)
    yz = mup + (yz * width)
    
    gauss_random_dl = cmplx(xz,yz,kind=8)
    
    return
    
  end function gauss_random_dl   

!*************************************************************************************************!

end module DL
