MODULE sb

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

  subroutine readparams_sb
  
    !   This subroutine reads the system specific parameters from the rundata.csv file
  
    implicit none
    character(LEN=100)::LINE1,LINE2,LINE3,LINE4,LINE5,LINE6,LINE7,LINE8,LINE9,LINE10,LINE11
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=128, file='rundata.csv', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening rundata.csv file file'
      errorflag = 1
      return
    end if

    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)
    read(128,*,iostat=ierr)LINE1,LINE2,LINE3,LINE4,LINE5,LINE6,LINE7,LINE8,LINE9,LINE10,LINE11
    
    close(128)

    read(LINE1,*,iostat=ierr)delta_sb
    if (ierr.ne.0) then
          write(0,"(a)") "Error reading Delta value"
          errorflag = 1
          return
        end if
        n = n+1
    read(LINE2,*,iostat=ierr)eps_sb
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading Epsilon value"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE3,*,iostat=ierr)beta_sb
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading beta value"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE4,*,iostat=ierr)wc_exp
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading wc value for exponential cutoff"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE5,*,iostat=ierr)kondo_exp
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading kondo parameter value"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE6,*,iostat=ierr)wmax_exp
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading wmax value for exponential cutoff"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE7,*,iostat=ierr)wc_dl
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading wc value for drude lorentz cutoff"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE8,*,iostat=ierr)lambda_dl
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading lambda value for drude lorentz cutoff"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE9,*,iostat=ierr)wmax_dl
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading wmax factor value for drude lorentz cutoff"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE10,*,iostat=ierr)uplimnorm
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading upper limit of the norm"
      errorflag = 1
      return
    end if
    n = n+1
    read(LINE11,*,iostat=ierr)lowlimnorm
    if (ierr.ne.0) then
      write(0,"(a)") "Error reading lower limit of the norm"
      errorflag = 1
      return
    end if
    n=n+1
    
    if (specden=="EXP") then 
      wmax_sb = wc_exp * wmax_exp
    else if (specden=="DL") then 
      wmax_sb = wc_dl * wmax_dl
    else if ((specden.ne."UBO").and.(specden.ne."LHC")) then
      write (0,"(a)") "Error in recognising spectral density while calculating wmax_sb"
      errorflag = 1
      return
    end if

    if (n.ne.11) then
      write(0,"(a)") "Not all required variables read in readparams_sb subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_sb

!--------------------------------------------------------------------------------------------------

  subroutine genwm_sb(wm_sb, Cm_sb, reps)   !   Level 2 Subroutine

    !   Generates the array of discretised frequencies over M degrees of freedom

    implicit none
    integer::m, ierr
    real(kind=8), dimension(:), allocatable, intent(inout) :: wm_sb, Cm_sb
    integer, intent(in) :: reps

    if (errorflag .ne. 0) return
    ierr = 0

    if (ndim.le.0) then
      write(0,"(a)") "Number of dimensions not correctly read or stored"
      errorflag = 1
      return
    else if (.not.allocated(wm_sb)) then
      write(0,"(a)") "genwm subroutine called but wm not allocated"
      errorflag=1
      return
    else if (.not.allocated(Cm_sb)) then
      write(0,"(a)") "genwm subroutine called but Cm not allocated"
      errorflag=1
      return
    else if (size(Cm_sb).ne.size(wm_sb)) then
      write(0,"(a)") "genwm subroutine called but Cm and wm do not have the same sizes"
      errorflag=1
      return
    else if (freqflg_sb.eq.0) then
      if (specden.eq."EXP") then
        do m=1,size(wm_sb)
          wm_sb(m)=-1.0d0*wc_exp*dlog(1.0d0-((dble(m)*(1.0d0-dexp(-1.0d0*wmax_sb/wc_exp)))/dble(ndim)))
          Cm_sb(m)=wm_sb(m)*sqrt(kondo_exp*wc_exp*(1.0d0-dexp(-1.0d0*wmax_sb/wc_exp))/dble(ndim))
        end do
      else if (specden.eq."DL") then
        do m=1,size(wm_sb)
          wm_sb(m)=wc_dl*dtan((dble(m)/dble(ndim))*datan(wmax_sb/wc_dl))
          Cm_sb(m)=wm_sb(m)*sqrt((4.0d0*lambda_dl/(sqrtpi*sqrtpi*wc_dl))*datan(wmax_sb/wc_dl)/dble(ndim))                           
        end do
      else
        write(0,"(3a)") "Error! Spectral Density set as ", specden, " but this is incompatible with frequency calculation."
        errorflag = 1
        return
      end if 
    else
      call readfreq(Cm_sb, wm_sb, reps)
    end if

    return

  end subroutine genwm_sb

!--------------------------------------------------------------------------------------------------

  subroutine gensig_sb(sig_sb, wm_sb)   !   Level 2 Subroutine

    ! Generates the distribution of states over the M degrees of freedom in the initial gaussian.

    implicit none
    real(kind=8), dimension (:), allocatable, intent (inout) :: sig_sb
    real(kind=8), dimension (:), intent (in) :: wm_sb
    integer::m, ierr

    if (errorflag .ne. 0) return
    ierr = 0

    if (ndim.le.0) then
      write(0,"(a)") 'Number of dimensions not correctly read or stored'
      errorflag = 1
      return
    else if (.not.allocated(sig_sb)) then
      write(0,"(a)") 'Allocation error in sig'
      errorflag = 1
      return
    else
      do m=1,size(sig_sb)
        sig_sb(m) = 1.0d0/(sqrt(dexp(beta_sb*wm_sb(m))-1.0d0))
      end do
    end if

    return

  end subroutine gensig_sb

!--------------------------------------------------------------------------------------------------

  subroutine genzinit_sb(mup, muq,reps)   !   Level 1 Subroutine

    ! Generates the initial Gaussian, used as a central point for the larger basis set (which is calculated in a different function)
    
    ! z_in is defined as being subject to the distribution rho(zin) = sigma_sb(m)*exp(-1.0*sigma_sb(m)*abs(zin(m)))

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq
    integer, intent(in) :: reps
    
    complex(kind=8),dimension(:,:), allocatable::H
    complex(kind=8), dimension (:), allocatable:: zin
    real(kind=8), dimension(:), allocatable :: sig_sb, wm_sb, Cm_sb
    real(kind=8) :: Ezin
    integer::m, n, recalc, ierr

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

    allocate (sig_sb(ndim), stat = ierr)
    if (ierr==0) allocate (wm_sb(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in sig or wm allocation in genzinit"
      errorflag=1
      return
    end if

    call genwm_sb(wm_sb, Cm_sb,reps)
    call gensig_sb(sig_sb, wm_sb)    

    do while (recalc == 1)
      do m=1,ndim
        if (randfunc.eq."ZBQL") then
          mup(m)=ZBQLNOR(mu,sig_sb(m)*sigp)
          muq(m)=ZBQLNOR(mu,sig_sb(m)*sigq)
          zin(m)=cmplx(muq(m),mup(m),kind=8)
        else
          zin(m)=gauss_random_sb(sig_sb(m),mu,mu)
          muq(m)=dble(zin(m))
          mup(m)=dimag(zin(m))
        end if
      end do
      if (ECheck.eq."YES") then
        call Hij_sb(H,zin,zin,wm_sb,Cm_sb)
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

    deallocate(zin, H, sig_sb, wm_sb, Cm_sb, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of zin, H, wm or sig in genzinit"
      errorflag=1
      return
    end if 

    return

  end subroutine genzinit_sb

!--------------------------------------------------------------------------------------------------
  
  subroutine Hord_sb(bs, H, t, reps)

    implicit none
    integer::k, j, r, s, ierr
    type(basisfn),dimension(:),intent(in)::bs
    type (hamiltonian), dimension (:,:), allocatable, intent(inout) :: H
    complex(kind=8), dimension (:,:), allocatable :: Hjk_mat
    real(kind=8), dimension(:), allocatable :: wm_sb, Cm_sb
    real(kind=8), intent (in) :: t
    integer, intent(in) :: reps

    if (errorflag .ne. 0) return

    allocate(Hjk_mat(npes,npes), stat = ierr)
    if (ierr==0) allocate (wm_sb(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of Hjk_mat, wm_sb and Cm_sb matrix in Hord_sb"
      errorflag=1
      return
    end if
       
    call genwm_sb(wm_sb, Cm_sb,reps)

    do k=1,size(H,2)
      do j=k,size(H,1)
        call Hij_sb(Hjk_mat,bs(j)%z,bs(k)%z, wm_sb, Cm_sb)
        do s=1,size(Hjk_mat,2)
          do r=1,size(Hjk_mat,1)
            H(j,k)%Hjk(r,s) = Hjk_mat(r,s)
            if (j.ne.k) then
              H(k,j)%Hjk(r,s) = dconjg(H(j,k)%Hjk(r,s))
            end if
          end do
        end do
      end do
    end do
    
    deallocate (Hjk_mat, wm_sb, Cm_sb, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of Hjk_mat, wm_sb or Cm_sb matrix in Hord_sb"
      errorflag=1
      return
    end if

    return

  end subroutine Hord_sb

!------------------------------------------------------------------------------------

  subroutine Hij_sb(H,z1,z2,wm_sb,Cm_sb)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    complex(kind=8):: Hbath, Hcoup
    real(kind=8), dimension(:), intent(in) :: wm_sb, Cm_sb
    real(kind=8) :: chk
    integer :: ierr

    if (errorflag .ne. 0) return

    ierr = 0

    Hbath = Hb_sb(z1,z2,wm_sb)
    Hcoup = Hc_sb(z1,z2,wm_sb, Cm_sb)

    chk = sum((abs(z1(1:ndim)-z2(1:ndim)))**2)

    if (chk.le.20000) then
      H(1,1) = Hbath+eps_sb+Hcoup
      H(1,2) = cmplx(delta_sb,0.0d0,kind=8)!+Hcoup
      H(2,1) = cmplx(delta_sb,0.0d0,kind=8)!+Hcoup
      H(2,2) = Hbath-eps_sb-Hcoup
    else
      H = (0.0d0, 0.0d0)
    end if

    return   

  end subroutine Hij_sb
  
!------------------------------------------------------------------------------------

  subroutine Hijdiag_sb(H,z,reps)

    implicit none
    complex(kind=8), dimension (:,:), intent(in)::z
    complex(kind=8), dimension(:,:,:), intent (inout)::H
    integer, intent(in) :: reps
    
    complex(kind=8), dimension(:), allocatable :: zin
    complex(kind=8):: Hbath, Hcoup
    real(kind=8), dimension(:), allocatable :: wm_sb, Cm_sb
    real(kind=8) :: chk
    integer :: k, m, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate (wm_sb(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
    if (ierr==0) allocate (zin(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in zin, wm or Cm allocation in Hijdiag"
      errorflag=1
      return
    end if

    call genwm_sb(wm_sb, Cm_sb,reps)
    
    do k=1,size(z,1)    
    
      do m=1,size(z,2)
        zin(m) = z(k,m)
      end do    

      Hbath = Hb_sb(zin,zin,wm_sb)
      Hcoup = Hc_sb(zin,zin,wm_sb, Cm_sb)
   
      H(k,1,1) = Hbath+eps_sb+Hcoup
      H(k,1,2) = cmplx(delta_sb,0.0d0,kind=8)!+Hcoup
      H(k,2,1) = cmplx(delta_sb,0.0d0,kind=8)!+Hcoup
      H(k,2,2) = Hbath-eps_sb-Hcoup

    end do

    deallocate(wm_sb, Cm_sb, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of wm or Cm in Hij"
      errorflag=1
      return
    end if 

    return   

  end subroutine Hijdiag_sb

!--------------------------------------------------------------------------------------------------

  function Hb_sb(z1,z2,wm_sb)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    real(kind=8), dimension(:), intent(in) :: wm_sb
    complex(kind=8)::Hb_sb
    integer :: m

    if (errorflag .ne. 0) return

    Hb_sb = (0.0d0, 0.0d0)

    do m=1,ndim
      Hb_sb = Hb_sb + (wm_sb(m)*(dconjg(z1(m))*z2(m)+0.5d0))
    end do

    return

  end function Hb_sb
    
!--------------------------------------------------------------------------------------------------

  function Hc_sb(z1,z2,wm_sb, Cm_sb)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    real(kind=8), dimension(:), intent(in) :: wm_sb, Cm_sb
    complex(kind=8)::Hc_sb
    integer :: m

    if (errorflag .ne. 0) return

    Hc_sb = (0.0d0, 0.0d0)

    do m=1,ndim
      Hc_sb = Hc_sb + ((Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))*(dconjg(z1(m))+z2(m)))
    end do

    return

  end function Hc_sb

!--------------------------------------------------------------------------------------------------

  function dh_dz_sb(z,reps)

    implicit none
    complex(kind=8),dimension(:,:),intent(in)::z
    integer, intent(in) :: reps
    
    complex(kind=8),dimension(size(z,1),npes,npes,size(z,2)) :: dh_dz_sb    
    real(kind=8), dimension(:), allocatable :: wm_sb, Cm_sb
    integer :: ierr, k

    if (errorflag .ne. 0) return

    if (npes.ne.2) then
      write(0,"(a)") "No. PES is not equal to 2, but Spin Boson 2 PES derivative called"
      errorflag = 1
      return
    end if

    allocate (wm_sb(ndim), stat=ierr)
    if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in wm or Cm allocation in dh_dz"
      errorflag=1
      return
    end if 

    call genwm_sb(wm_sb, Cm_sb, reps)

    do k=1,size(z,1)
      dh_dz_sb(k,1,1,:) = dhdz_sb_11(z(k,:),wm_sb,Cm_sb)
      dh_dz_sb(k,1,2,:) = dhdz_sb_12(wm_sb,Cm_sb)
      dh_dz_sb(k,2,1,:) = dhdz_sb_21(wm_sb,Cm_sb)
      dh_dz_sb(k,2,2,:) = dhdz_sb_22(z(k,:),wm_sb,Cm_sb)
    end do

!    do k=1,size(z,1)
!      dh_dz_sb(1,1,:) = dhdz_sb_11(z(k,:),wm_sb)
!      dh_dz_sb(1,2,:) = dhdz_sb_12(wm_sb,Cm_sb)
!      dh_dz_sb(2,1,:) = dhdz_sb_21(wm_sb,Cm_sb)
!      dh_dz_sb(2,2,:) = dhdz_sb_22(z(k,:),wm_sb)
!    end do

    deallocate (wm_sb, Cm_sb, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in wm or Cm allocation in dh_dz"
      errorflag=1
      return
    end if 

    return

  end function dh_dz_sb

!--------------------------------------------------------------------------------------------------

  function dhdz_sb_11(z,wm_sb, Cm_sb)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_sb_11
    complex(kind=8),dimension(:),intent(in)::z
    real(kind=8),dimension(:),intent(in)::wm_sb, Cm_sb
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
!      dhdz_sb_11(m) = wm_sb(m)*z(m)
      dhdz_sb_11(m) = wm_sb(m)*z(m) + (Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))
    end do

    return

  end function dhdz_sb_11

!--------------------------------------------------------------------------------------------------

  function dhdz_sb_12(wm_sb,Cm_sb)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_sb_12
    real(kind=8),dimension(:),intent(in)::wm_sb, Cm_sb
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
!      dhdz_sb_12(m) = (Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))
      dhdz_sb_12(m) = (0.0d0,0.0d0)
    end do

    return

  end function dhdz_sb_12

!--------------------------------------------------------------------------------------------------

  function dhdz_sb_21(wm_sb,Cm_sb)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_sb_21
    real(kind=8),dimension(:),intent(in)::wm_sb, Cm_sb
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
!      dhdz_sb_21(m) = (Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))
      dhdz_sb_21(m) = (0.0d0,0.0d0)
    end do

    return

  end function dhdz_sb_21

!--------------------------------------------------------------------------------------------------

  function dhdz_sb_22(z,wm_sb,Cm_sb)

    implicit none
    complex(kind=8),dimension(ndim)::dhdz_sb_22
    complex(kind=8),dimension(:),intent(in)::z
    real(kind=8),dimension(:),intent(in)::wm_sb, Cm_sb
    integer::m

    if (errorflag .ne. 0) return

    do m=1,ndim
      dhdz_sb_22(m) = wm_sb(m)*z(m) - (Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))
!      dhdz_sb_22(m) = wm_sb(m)*z(m)
    end do

    return

  end function dhdz_sb_22

!--------------------------------------------------------------------------------------------------

  function gauss_random_sb (width, muq, mup)
  
    implicit none
    complex(kind=8) :: gauss_random_sb
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
    
    gauss_random_sb = cmplx(xz,yz,kind=8)
    
    return
    
  end function gauss_random_sb 
  
!--------------------------------------------------------------------------------------------------

  subroutine readfreq(Cm_sb, wm_sb, reps) 
  
    implicit none
    
    real(kind=8), dimension(:), allocatable, intent(inout) :: wm_sb, Cm_sb
    integer, intent(in) :: reps
    
    real :: r
    integer(kind=8)::ranseed
    integer,allocatable::ranseed_array(:)
    integer:: ranseed_size, rint
    integer::i,clock, m, ierr
    character(LEN=100)::LINE
    character(LEN=3)::repstr
    
    if (errorflag .ne. 0) return
    
    ierr=0
    
    call random_seed(size=ranseed_size)
    allocate(ranseed_array(ranseed_size))
    call system_clock(count=clock)
    ranseed_array = clock + 37* (/ (i-1,i=1,ranseed_size) /)
    call random_seed(put=ranseed_array)
    deallocate(ranseed_array)
    CALL RANDOM_NUMBER(r)
    r=r*3232768.0*dble(reps)
    rint = int(r)
    
    write(repstr,"(i3.3)") reps 
      
    open(unit=rint, file="freq"//trim(repstr)//".dat", status='old', iostat=ierr)
    
    if (ierr.ne.0) then
      write(0,"(3a)") 'Error in opening freq'//trim(repstr)//'.dat file'
      errorflag = 1
      return
    end if
    
    do i=1,size(wm_sb)
     
      read(rint,*,iostat=ierr) m, wm_sb(m), Cm_sb(m)
      if (i.ne.m) then
        write (6,"(a,i0,a,i0)") "Mismatch in frequency reading. Read entry ",m," but expected ",i
        errorflag = 1
        return
      end if
                 
    end do
    
    close(rint)
      
    return    
    
  end subroutine readfreq   

!*************************************************************************************************!

end module sb
