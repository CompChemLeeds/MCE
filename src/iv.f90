MODULE iv

  use globvars

!***********************************************************************************!
!*
!*    Inverted Gaussian Module
!*           
!*    Contains subroutines for:
!*
!*    1) Reading the Morse Potential Parameters
!*    2) Calculating the real and imaginary parts of the m dimensional initial 
!*        wavefunction zinit
!*    3) Calculating a single element of the Hamiltonian matrix
!*    4) Calculating the derivative of the hamiltonian
!*      
!***********************************************************************************!

contains

!------------------------------------------------------------------------------------

  subroutine readparams_iv
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
      if(LINE=='IVm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,mass_iv
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading mass value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='IVw') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,freq_iv
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading frequency value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='IVlambda') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,lambda_iv
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading lambda value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='IVIntensity') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,inten_iv
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading laser intensity value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='IVupnorm') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,uplimnorm
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading upper limit of the norm"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='IVdownnorm') then
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
      write(0,"(a)") "Not all required variables read in readparams_iv subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readparams_iv

!------------------------------------------------------------------------------------

  subroutine genzinit_iv(mup, muq)   !   Level 1 Subroutine

    implicit none
    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

    if (errorflag .ne. 0) return

    muq(1:ndim) = 0.0d0*sigp
    mup(1:ndim) = 0.0d0*sigq

    return

  end subroutine genzinit_iv

!------------------------------------------------------------------------------------
  
  subroutine Hord_iv(bs, H, t)

    implicit none
    integer::k, j, r, s, ierr
    type(basisfn),dimension(:),intent(in)::bs
    type (hamiltonian), dimension (:,:), allocatable, intent(inout) :: H
    complex(kind=8), dimension (:,:), allocatable :: Hjk_mat
    real(kind=8), intent (in) :: t

    if (errorflag .ne. 0) return

    allocate(Hjk_mat(npes,npes), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of Hjk_mat matrix in Hord"
      errorflag=1
      return
    end if

    do k=1,size(H,2)
      do j=k,size(H,1)
        call Hij_iv(Hjk_mat,bs(j)%z,bs(k)%z,t)
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
    
    deallocate (Hjk_mat, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of Hjk_mat matrix in Hord"
      errorflag=1
      return
    end if

    return

  end subroutine Hord_iv

!------------------------------------------------------------------------------------

  subroutine Hij_iv(H,z1,z2,t)

    implicit none 
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8), dimension(:,:), intent (inout)::H
    real(kind=8), intent (in) :: t
    integer :: m, ierr
    complex(kind=8), dimension (:), allocatable :: Htemp, z1c, rho
    real (kind=8) :: eta, field

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Inverse Gaussian"
      errorflag = 1
      return
    end if

    allocate (Htemp(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating Htemp in Hij_iv"
      errorflag=1
    end if

    eta = lambda_iv*gam/(gam+lambda_iv)

    allocate(z1c(size(z1)))
    allocate(rho(size(z1)))
    z1c(1:ndim)=dconjg(z1(1:ndim))

    do m=1,ndim
      rho(m) = (z1c(m)+z2(m))/sqrt(2.0*gam)
    end do
    field = inten_iv*dcos(freq_iv*t)
!    field = inten_iv*dsin(freq_iv*t)*(dsin(freq_iv*t/4.0d0))**2.0d0 

    do m=1,ndim
      Htemp(m) = (0.0d0,0.0d0)
      Htemp(m) = Htemp(m) - ((hbar**2)*gam/(4.0d0*mass_iv))&
                   *(z1c(m)**2.0d0+z2(m)**2.0d0-2.0d0*z1c(m)*z2(m)-1.0d0)
      if (m==1) Htemp(m) = Htemp(m) + rho(m)*field
    end do

    H(1,1) = sum(Htemp(1:ndim))
    
    do m=1,ndim
      Htemp(m) = (0.0d0,0.0d0)
      Htemp(m) = Htemp(m) - sqrt(eta/lambda_iv)*cdexp(-1.0d0*eta*rho(m)**2.0d0)
    end do
    
    H(1,1) =  H(1,1) + product(Htemp(1:ndim))  
    
    if (H(1,1)/=H(1,1)) then
      write(0,"(a)") "Error! Hamiltonian element NaN"
      errorflag = 1
      return
    end if

    deallocate(Htemp, z1c)

    return   

  end subroutine Hij_iv
  
!------------------------------------------------------------------------------------

  subroutine Hijdiag_iv(H,z,t)

    implicit none 
    complex(kind=8), dimension (:,:), intent(in)::z
    complex(kind=8), dimension(:,:,:), intent (inout)::H
    real(kind=8), intent (in) :: t
    integer :: k, m, ierr
    complex(kind=8), dimension (:), allocatable :: Htemp, rho
    complex(kind=8) :: zc
    real (kind=8) :: eta, field

    if (errorflag .ne. 0) return

    if (npes.ne.1) then
      write(0,"(a)") "Error! There is more than 1 pes for the Inverse Gaussian"
      errorflag = 1
      return
    end if

    allocate (Htemp(ndim), stat=ierr)
    if (ierr==0) allocate(rho(size(z)), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating Htemp and rho in Hij_iv"
      errorflag=1
    end if

    eta = lambda_iv*gam/(gam+lambda_iv)
    field = inten_iv*dcos(freq_iv*t)
!    field = inten_iv*dsin(freq_iv*t)*(dsin(freq_iv*t/4.0d0))**2.0d0     

    do k=1,size(H,1)

      do m=1,ndim
        zc=dconjg(z(k,m))
        rho(m) = (zc+z(k,m))/sqrt(2.0*gam)
        Htemp(m) = (0.0d0,0.0d0)
        Htemp(m) = Htemp(m) - ((hbar**2)*gam/(4.0d0*mass_iv))&
                     *(zc**2.0d0+z(k,m)**2.0d0-2.0d0*zc*z(k,m)-1.0d0)
        if (m==1) Htemp(m) = Htemp(m) + rho(m)*field
      end do
  
      H(k,1,1) = sum(Htemp(1:ndim))
      
      do m=1,ndim
        Htemp(m) = (0.0d0,0.0d0)
        Htemp(m) = Htemp(m) - sqrt(eta/lambda_iv)*cdexp(-1.0d0*eta*rho(m)**2.0d0)
      end do
      
      H(k,1,1) =  H(k,1,1) + product(Htemp(1:ndim))  
    
      if (H(k,1,1)/=H(k,1,1)) then
        write(0,"(a,i0,a)") "Error! Hamiltonian element ", k, " NaN"
        errorflag = 1
        return
      end if

    end do

    deallocate(Htemp, stat=ierr)
    if (ierr==0) deallocate(rho, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating Htemp and rho in Hij_iv"
      errorflag=1
    end if

    return   

  end subroutine Hijdiag_iv

!------------------------------------------------------------------------------------

  function dh_dz_iv(z,t)

    implicit none
    complex(kind=8),dimension(:,:),intent(in)::z
    real(kind=8), intent (in) :: t
    
    complex(kind=8),dimension(size(z,1),npes,npes,size(z,2)) :: dh_dz_iv
    complex(kind=8), dimension(:), allocatable :: zc, rho  
    complex(kind=8) :: dhdztmp, rho2
    real(kind=8) :: eta, field
    integer :: m, k, ierr

    if (errorflag .ne. 0) return

    allocate (zc(size(z,2)), stat=ierr)
    if (ierr==0) allocate (rho(size(z,2)), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating zc and rho arrays in dh_dz_iv function"
      errorflag = 1
      return
    end if
    
    eta = lambda_iv*gam/(gam+lambda_iv)

    field = inten_iv*dcos(freq_iv*t)
!    field = inten_iv*dsin(freq_iv*t)*(dsin(freq_iv*t/4.0d0))**2.0d0    
    
    do k=1,size(z,1)
      do m=1,size(z,2)
        zc(m)=dconjg(z(k,m))
        rho(m) = (zc(m)+z(k,m))/sqrt(2.0*gam)
      end do
      
      rho2 = sum(rho(:)**2.0d0)

      do m=1,size(z,2)   
        dhdztmp = (0.0d0,0.0d0) 
        dhdztmp = dhdztmp - ((hbar**2)*gam/(2.0d0*mass_iv))*(zc(m)-z(k,m))
        dhdztmp = dhdztmp + eta*sqrt((2.0/gam)*((eta/lambda_iv)**dble(ndim)))*rho(m)*&
                            exp(-1.0d0*eta*rho2)
        if (m==1) dhdztmp = dhdztmp + (1.0/sqrt(2.0d0*gam))*field
        dh_dz_iv (k,1,1,m) = dhdztmp
      end do
    end do
    
    deallocate (zc, stat=ierr)
    if (ierr==0) deallocate (rho, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating zc and rho arrays in dh_dz_iv function"
      errorflag = 1
      return
    end if    

    return

  end function dh_dz_iv

!------------------------------------------------------------------------------------

  function dipole_iv(bs)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8) :: dipole_iv, ovrlp, rho2
    complex(kind=8), dimension(:), allocatable :: D, Dc, zk, zj, rho
    real(kind=8), dimension (:), allocatable :: s
    real(kind=8) :: eta
    integer::k,j,m,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    dipole_iv = (0.0d0,0.0d0)
    eta = gam*lambda_iv/(gam+lambda_iv)

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
        rho2 = (0.0d0,0.0d0)
        do m=1,ndim
          zk(m) = bs(k)%z(m)
          zj(m) = bs(j)%z(m)
          rho(m) = (dconjg(zj(m))+zk(m))/sqrt(2.0*gam)
        end do
        rho2 = sum(rho(1:ndim)**2.0d0)
        ovrlp = product(cdexp((dconjg(zj(1:ndim))*zk(1:ndim))&
              -(0.5d0*dconjg(zj(1:ndim))*zj(1:ndim))&
              -(0.5d0*dconjg(zk(1:ndim))*zk(1:ndim))))
        dipole_iv = dipole_iv + (-1.0d0*Dc(j)*D(k)*cdexp(i*(s(k)-s(j)))* ovrlp*&
                                    2.0d0*eta*sqrt((eta/lambda_iv)**dble(ndim))*rho(1)&
                                                *cdexp(-1.0d0*eta*rho2))
      end do
    end do

    deallocate (D, stat = ierr)
    if (ierr==0) deallocate (Dc, stat=ierr)
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

    return

  end function dipole_iv

!***********************************************************************************!

end module iv
