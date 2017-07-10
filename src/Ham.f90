MODULE Ham

  use globvars
  use alarrays
  use redirect

!***********************************************************************************!
!*
!*         Generalised Hamiltonian Module
!*           
!*   Contains subroutines for:
!*
!*      1) Calculating the overlap matrix with elements <z_j|z_k>
!*      2) Calculating the individual elements of the matrix calculated in 1)
!*      3) Calculating the overlap matrix with elements <phi_j|phi_k>
!*      4) Calculating the norm of the overall wavefunction
!*      5) Calculating the sum of the norms of the individual configurations
!*      6) Calculating the population in a particular electronic state
!*      7) Calculating the auto-correlation function of the wavefunction
!*      8) Calculating the Hamiltonian matrix
!*      9) Calculating the Ehrenfest Hamiltonian for an individual basis function
!*     10) Linear equation calculation used in calculating intial D or dD/dt.
!*     11) Matrix Inversion subroutine to be used as an alternative to 10)
!*         but increases runtime by ~3000% so not recommended for use.
!*      
!***********************************************************************************!

  contains

!***********************************************************************************!
!             Overlap Functions
!***********************************************************************************!

  function ovrlpmat(bs)

    implicit none
    type(basisfn), dimension (:), intent(in)::bs
    complex(kind=8), dimension(size(bs),size(bs))::ovrlpmat
    integer::k,j,n

    if (errorflag .ne. 0) return

    ovrlpmat = (0.0d0, 0.0d0)

    n = size(bs)

    do k=1,n
      do j=k,n
        if (j==k) then
          ovrlpmat(j,k)=(1.0d0,0.0d0)
        else
          ovrlpmat(j,k)=ovrlpij(bs(j)%z,bs(k)%z)
          ovrlpmat(k,j)=dconjg(ovrlpmat(j,k))
        end if
      end do
    end do
    
    return

  end function ovrlpmat

!------------------------------------------------------------------------------------

  function ovrlpij(z1,z2)

    implicit none
    complex(kind=8), dimension (:), intent(in)::z1,z2
    complex(kind=8)::ovrlpij
    real(kind=8)::chk
    integer :: m

    if (errorflag .ne. 0) return

    chk = sum((abs(z1(1:ndim)-z2(1:ndim)))**2)

    ovrlpij = (0.0d0,0.0d0)

    if (chk.le.20000) then
      do m=1,ndim
        ovrlpij  = ovrlpij  + ((dconjg(z1(m))*z2(m))-(0.5d0*dconjg(z1(m))*z1(m))&
                  -(0.5d0*dconjg(z2(m))*z2(m)))
      end do
      ovrlpij  = cdexp(ovrlpij )
    end if

    return

  end function ovrlpij

!------------------------------------------------------------------------------------

  function ovrlpphimat(bs)

    implicit none
    type(basisfn), dimension (:), intent(in)::bs
    complex(kind=8), dimension(size(bs),size(bs))::ovrlpphimat
    complex(kind=8)::asum
    integer::k,j,r,n

    if (errorflag .ne. 0) return

    n = size(bs)

    do k=1,n
      do j=k,n
          asum = (0.0d0,0.0d0)
          do r=1,npes
            asum = asum + (dconjg(bs(j)%a_pes(r))*bs(k)%a_pes(r))
          end do
          ovrlpphimat(j,k)=ovrlpij(bs(j)%z,bs(k)%z)*asum
        if (k.ne.j) then
          ovrlpphimat(k,j)=dconjg(ovrlpphimat(j,k))
        end if
      end do
    end do
    
  end function ovrlpphimat

!***********************************************************************************!
!             Norm and Population Functions
!***********************************************************************************!

  function norm(bs,ovrlp)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8), dimension(:,:), intent(in) :: ovrlp    
    complex(kind=8)::norm, asum
    real(kind=8)::absnorm
    integer::k, j, r, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    norm = (0.0d0, 0.0d0)
   
    do k=1,size(bs)
      do j=1,size(bs)
        asum = (0.0d0, 0.0d0)
        do r=1,npes
          asum = asum + (dconjg(bs(j)%a_pes(r))*bs(k)%a_pes(r))
        end do
        norm = norm + dconjg(bs(j)%D_big)*ovrlp(j,k)*asum*bs(k)%D_big
      end do
    end do

    absnorm=sqrt(dble(norm*dconjg(norm)))

    if ((absnorm.gt.1.2d0).and.(method.ne."AIMC1")) then
      write(0,"(a,a)") "Norm is greater than 1.2. Simulation has failed and",&
                        " trajectories are likely to explode"
      write(0,"(a,e16.8e3)") "Norm had a value of ", absnorm
      errorflag = 1
    end if

    return

  end function norm

!--------------------------------------------------------------------------------------------------

  function norm2(bs)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8)::norm2
    complex(kind=8), dimension(:), allocatable :: normar
    integer::r,j,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(normar(size(bs)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of normar array in norm2"
      errorflag=1
      return
    end if   

    norm2 = (0.0d0, 0.0d0)
    normar = (0.0d0, 0.0d0)

    if (method/="MCEv2") then
      write(0,"(a)") "Method is not valid for calculation of single configuration norm."
      write(0,"(a)") "This subroutine has been incorrectly called"
      errorflag=1
      return
    end if
    
    do j=1,size(bs)
      do r=1,npes
        normar(j) = normar(j) + dconjg(bs(j)%d_pes(r))*bs(j)%d_pes(r)
      end do
    end do

    norm2 = product(normar(1:size(bs)))

    deallocate(normar, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of normar array in norm2"
      errorflag=1
      return
    end if

    return

  end function norm2

!------------------------------------------------------------------------------------

  function pop(bs, pes, ovrlp)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8), dimension(:,:), intent(in) :: ovrlp
    complex(kind=8)::cpop
    integer::k,j,ierr
    integer, intent(in):: pes
    real(kind=8)::pop

    if (errorflag .ne. 0) return

    ierr = 0

    cpop = (0.0d0, 0.0d0)
  
    do k=1,size(bs)
      do j=1,size(bs)
        cpop = cpop + dconjg(bs(j)%D_big) * ovrlp(j,k) * bs(k)%D_big &
                * dconjg(bs(j)%d_pes(pes)) * bs(k)%d_pes(pes) &
                * cdexp(i*(bs(k)%s_pes(pes)-bs(j)%s_pes(pes)))
      end do
    end do

    pop = abs(cpop) 

    return

  end function pop

!------------------------------------------------------------------------------------

  function innorm(bs)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs    
    complex(kind=8)::innorm, asum
    complex(kind=8), dimension(:), allocatable::zj,zk 
    real(kind=8)::absnorm
    integer::k, j, r, m, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    innorm = (0.0d0, 0.0d0)

    allocate(zj(ndim), zk(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating z arrays"
      errorflag = 1
      return
    end if
    
    if (((basis=="GRID").or.(basis=="GRSWM")).and.(mod(in_nbf,2)==1)) then
      innorm = (1.0d0,0.0d0)
    else
      do k=1,size(bs)
        do j=1,size(bs)
          asum = (0.0d0, 0.0d0)
          do m=1,ndim
            zj(m)=bs(j)%z(m)
            zk(m)=bs(k)%z(m)
          end do
          do r=1,npes
            asum = asum + (dconjg(bs(j)%a_pes(r))*bs(k)%a_pes(r))
          end do
            innorm = innorm + dconjg(bs(j)%D_big)*ovrlpij(zj,zk)*asum*bs(k)%D_big
        end do
      end do
    end if

    absnorm=sqrt(dble(innorm*dconjg(innorm)))

!    if ((absnorm.gt.1.2d0).and.(cmprss=="N")) then
!      write(0,"(a,a)") "Norm is greater than 1.2. Simulation has failed and trajectories ",&
!                "are likely to explode"
!      errorflag = 1
!    end if

    deallocate(zj, zk, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating z arrays"
      errorflag = 1
      return
    end if

    return

  end function innorm

!------------------------------------------------------------------------------------

  function inpop(bs, pes)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8), dimension(:), allocatable::zj,zk 
    complex(kind=8)::cpop
    integer::k,j,ierr, m
    integer, intent(in):: pes
    real(kind=8)::inpop

    if (errorflag .ne. 0) return

    ierr = 0

    cpop = (0.0d0, 0.0d0)

    allocate(zj(ndim), zk(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error allocating z arrays"
      errorflag = 1
      return
    end if

    if (((basis=="GRID").or.(basis=="GRSWM")).and.(mod(in_nbf,2)==1)) then
      if (pes==in_pes) then
        inpop = 1.0d0
      else
        inpop = 0.0d0
      end if
    else    
      do k=1,size(bs)
        do j=1,size(bs)
          do m=1,ndim
            zj(m)=bs(j)%z(m)
            zk(m)=bs(k)%z(m)
          end do
          cpop = cpop + dconjg(bs(j)%D_big) * ovrlpij(zj,zk) * bs(k)%D_big &
                  * dconjg(bs(j)%d_pes(pes)) * bs(k)%d_pes(pes) * &
                    cdexp(i*(bs(k)%s_pes(pes)-bs(j)%s_pes(pes)))
        end do
      end do
      inpop = abs(cpop) 
    end if

    deallocate(zj, zk, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating z arrays"
      errorflag = 1
      return
    end if

    return

  end function inpop

!------------------------------------------------------------------------------------

  function acf(bs, mup, muq)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    real(kind=8), dimension(:), intent(in)::mup,muq
    real(kind=8), dimension(:), allocatable :: s
    complex(kind=8) :: acf
    complex(kind=8), dimension (:), allocatable :: zinit, Dbig, d, ovrlp
    integer::k,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(zinit(ndim), stat=ierr)
    if (ierr==0) allocate (Dbig(size(bs)), stat=ierr)
    if (ierr==0) allocate (d(size(bs)), stat=ierr)
    if (ierr==0) allocate (s(size(bs)), stat=ierr)
    if (ierr==0) allocate (ovrlp(size(bs)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of arrays in acf"
      errorflag=1
      return
    end if

    acf = (0.0d0,0.0d0)

    zinit(1:ndim) = cmplx(muq(1:ndim), mup(1:ndim), kind=8)

    do k=1,size(bs)
      Dbig(k)=bs(k)%D_big
      d(k)=bs(k)%d_pes(in_pes)
      s(k)=bs(k)%s_pes(in_pes)
      ovrlp(k)=ovrlpij(zinit, bs(k)%z)
    end do 

    do k=1,size(bs)
      acf = acf + (ovrlp(k)*Dbig(k)*d(k)*cdexp(i*s(k)))
    end do

    deallocate(zinit, Dbig, d, s, ovrlp, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of arrays in acf"
      errorflag=1
      return
    end if

    return

  end function acf

!------------------------------------------------------------------------------------

  function acfdim(bs, mup, muq)   !   Level 1 Function

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    real(kind=8), dimension(:), intent(in)::mup,muq
    real(kind=8), dimension(:), allocatable :: s
    real(kind=8) , dimension(3*ndim+3) :: acfdim
    complex(kind=8), dimension (:), allocatable :: zinit, Dbig, d, acft, ovrlpprod
    complex(kind=8), dimension (:,:), allocatable :: ovrlp
    real(kind=8) :: chk
    integer::k,m,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(zinit(ndim), stat = ierr)
    if (ierr==0) allocate(Dbig(size(bs)), stat=ierr)
    if (ierr==0) allocate(d(size(bs)), stat=ierr)
    if (ierr==0) allocate(s(size(bs)), stat=ierr)
    if (ierr==0) allocate(ovrlp(size(bs),ndim), stat=ierr)
    if (ierr==0) allocate(acft(ndim+1), stat=ierr)
    if (ierr==0) allocate(ovrlpprod(size(bs)), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of arrays in acf"
      errorflag=1
      return
    end if

    acft = (0.0d0,0.0d0)
    ovrlpprod = (1.0d0,0.0d0)

    zinit(1:ndim) = cmplx(muq(1:ndim), mup(1:ndim), kind=8)

    do k=1,size(bs)
      Dbig(k)=bs(k)%D_big
      d(k)=bs(k)%d_pes(in_pes)
      s(k)=bs(k)%s_pes(in_pes)
      do m=1,ndim
        chk = (abs(zinit(m)-bs(k)%z(m)))**2
        if (chk.gt.20000d0) then
          ovrlp(k,m) = (0.0d0,0.0d0)
        else
          ovrlp(k,m) = cdexp((dconjg(zinit(m))*bs(k)%z(m))-&
             (0.5d0*dconjg(zinit(m))*zinit(m))-(0.5d0*dconjg(bs(k)%z(m))*bs(k)%z(m)))
          ovrlpprod(k) = ovrlpprod(k) * ovrlp(k,m) 
        end if       
      end do
    end do 

    do m=1,ndim
      do k=1,size(bs)
        acft(m) = acft(m) + (ovrlp(k,m)*Dbig(k)*d(k)*cdexp(i*s(k)))
      end do
    end do

    do k=1,size(bs)
      acft(ndim+1) = acft(ndim+1) + (ovrlpprod(k)*Dbig(k)*d(k)*cdexp(i*s(k)))
    end do

    do m=1,ndim
      acfdim(m+1) = dimag(i*acft(m))
      acfdim(m+ndim+2) = dimag(acft(m))
      acfdim(m+2*ndim+3) = abs(acft(m))
    end do

    acfdim(1) = dimag(i*acft(ndim+1))
    acfdim(ndim+2) = dimag(acft(ndim+1))
    acfdim(2*ndim+3) = abs(acft(ndim+1))

    deallocate(zinit, Dbig, d, s, ovrlp, acft, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of arrays in acf"
      errorflag=1
      return
    end if

    return

  end function acfdim


!***********************************************************************************!
!            Hamiltonian Subroutines and Functions
!***********************************************************************************!

  subroutine Hord(bs, H, t)

    implicit none
    integer::k, j, r, s, ierr
    type(basisfn),dimension(:),intent(in)::bs
    type (hamiltonian), dimension (:,:), allocatable, intent(inout) :: H
    complex(kind=8), dimension (:,:), allocatable :: Hjk_mat
    real(kind=8), intent (in) :: t

    if (errorflag .ne. 0) return

    if (.not.allocated(H)) then
      write(0,"(a)") "Error! Hord has not been allocated in Hord subroutine."    
      errorflag = 1
      return
    end if

    allocate(Hjk_mat(npes,npes), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of Hjk_mat matrix in Hord"
      errorflag=1
      return
    end if

    do k=1,size(H,2)
      do j=k,size(H,1)
        call Hij(Hjk_mat,bs(j)%z,bs(k)%z,t)
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

  end subroutine Hord

!------------------------------------------------------------------------------------

  function Hehr(bf,t)   !   Level 1 Function

    implicit none
    type (basisfn), intent(inout)::bf
    real(kind=8), intent (in) :: t
    complex(kind=8)::Hehr, asum
    complex(kind=8),dimension(:,:), allocatable::H
    integer::r,s,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(H(npes,npes), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of H matrix in Hehr"
      errorflag=1
      return
    end if

    Hehr = (0.0d0, 0.0d0)
    asum = (0.0d0, 0.0d0)

    call Hij(H,bf%z,bf%z,t) 

    do r=1,npes
      do s=1,npes
        Hehr = Hehr + H(r,s)*dconjg(bf%a_pes(r))*bf%a_pes(s)
      end do
      asum = asum + dconjg(bf%a_pes(r))*bf%a_pes(r)
    end do

    Hehr = Hehr/asum

    deallocate(H, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of H matrix in Hehr"
      errorflag=1
      return
    end if

    return

  end function Hehr

!***********************************************************************************!
!              Subroutines for Calculation of matrix inverse
!
!      -lineq calculates the solution to the system of linear equations Ovrlp*D=C
!      -matinv calculates the inverse of the matrix from eigenvalues and eigenvectors
!                 and then solves Ovrlp*D=C through D=Ovrlp^(-1)*C
!  
!***********************************************************************************!

  subroutine lineq(ovrlp, C_k, D)   !   Level 2 Subroutine

    implicit none
    complex(kind=8),dimension(:),intent(inout)::D, C_k
    complex(kind=8),dimension(:,:),intent(inout):: ovrlp
    integer::nrhs ,info, n, ierr
    integer, dimension(:), allocatable::IPIV

    if (errorflag .ne. 0) return
    ierr = 0
    nrhs = 1
    info = 0

    allocate(IPIV(size(C_k)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error IPIV allocation for Linear Algebra function"
      errorflag=1
      return
    end if

    n = size(C_k)
    if ((size(D).ne.n).or.(size(ovrlp,1).ne.n).or.(size(ovrlp,2).ne.n)) then
      write(0,"(a)") "Size mismatch in arrays for zgesv subroutine"
    end if

    call zgesv(n,nrhs,ovrlp,n,IPIV,C_k,n,info)

    if (info.ne.0) then
      write(0,"(a,i0)")"Error in zgesv. Info value is ", info
      errorflag=1
      return
    end if

    D = C_k

    deallocate(IPIV, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of local arrays in lineq"
      errorflag=1
      return
    end if

    return

  end subroutine lineq

!------------------------------------------------------------------------------------

  subroutine matinv(ovrlp, C_k, D, restart)   !   Level 2 Subroutine

    implicit none
    complex(kind=8),dimension(:),intent(inout)::D, C_k
    complex(kind=8),dimension(:,:),intent(inout):: ovrlp
    complex(kind=8),dimension(:,:),allocatable::om_inv
    complex(kind=8),dimension(:),allocatable::WORK
    real(kind=8),dimension(:),allocatable::RWORK, om_eigen
    integer, intent(inout) :: restart
    integer::l, j, k, n, info, LWORK, nonherm, ierr

    if (errorflag .ne. 0) return
    ierr = 0
    info = 0
    nonherm = 0
    n = size(C_k)


    LWORK = 2*n
    allocate(om_eigen(n), stat = ierr)
    allocate(WORK(LWORK), stat = ierr)
    allocate(RWORK(3*n-2), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of arrays for matrix inversion"
      errorflag=1
      return
    end if

    if ((size(D).ne.n).or.(size(ovrlp,1).ne.n).or.(size(ovrlp,2).ne.n)) then
      write(0,"(a)") "Error in matinv subroutine! Array sizes are mismatched."
      errorflag = 1
      return
    end if

    call zheev('V','U',n,ovrlp,n,om_eigen,WORK,LWORK,RWORK,info)

    if (info.ne.0) then
      write(0,"(a,i0)")"Error in zheev. Info value is ", info
      errorflag=1
      return
    end if

    if (minval(abs(om_eigen)) == 0.0d0) then
      write(0,"(a)") "Error in inverse matrix operation. There exists a zero eigenvalue"
      restart = 1
      return       
    end if

    allocate(om_inv(n,n), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in inverse matrix allocation"
      errorflag=1
      return
    end if
    om_inv=(0.0d0,0.0d0)

    do l=1,n
      do j=1,n
        do k=1,n
          if (abs(om_eigen(k))/maxval(abs(om_eigen)).gt.1.0d-3) then
            om_inv(l,j)=om_inv(l,j) + ovrlp(l,k)*dconjg(ovrlp(j,k))/om_eigen(k)
          end if
        end do
      end do
    end do

    do k=1,n
      do j=1,k-1
        if (dconjg(om_inv(k,j)).ne.om_inv(j,k)) then
          nonherm = nonherm + 1
        end if
      end do
    end do

    if (nonherm.ne.0) then
      write(0,"(a,i5,a,i5,a,i5)")"Error! Inverse matrix is Non-Hermitian. ",2*nonherm,&
                " unmatched off-diagonal elements found in a ", n,"x",n," matrix"
      restart = 1
      return
    end if

    D = matmul(om_inv,C_k)
    deallocate(om_inv, stat = ierr)
    deallocate(om_eigen, stat = ierr)
    deallocate(WORK, stat = ierr)
    deallocate(RWORK, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of local arrays in matinv"
      errorflag=1
      return
    end if

    return

  end subroutine matinv

!------------------------------------------------------------------------------------

  subroutine matinv2(ovrlp, C_k, D)

    implicit none
    complex(kind=8),dimension(:),intent(inout)::D, C_k
    complex(kind=8),dimension(:,:),intent(inout):: ovrlp
    complex(kind=8),dimension(:,:),allocatable::om_inv
    complex(kind=8),dimension(:),allocatable::WORK
    real(kind=8),dimension(:),allocatable::RWORK, om_eigen
    integer::l, j, k, n, info, LWORK, nonherm, ierr

    if (errorflag .ne. 0) return
    ierr = 0
    info = 0
    nonherm = 0
    n = size(C_k)

    LWORK = 2*n
    allocate(om_eigen(n), stat = ierr)
    allocate(WORK(LWORK), stat = ierr)
    allocate(RWORK(3*n-2), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of arrays for matrix inversion"
      errorflag=1
      return
    end if

    if ((size(D).ne.n).or.(size(ovrlp,1).ne.n).or.(size(ovrlp,2).ne.n)) then
      write(0,"(a)") "Error in matinv subroutine! Array sizes are mismatched."
      errorflag = 1
      return
    end if

    call zheev('V','U',n,ovrlp,n,om_eigen,WORK,LWORK,RWORK,info)

    if (info.ne.0) then
      write(0,"(a,i0)")"Error in zheev. Info value is ", info
      errorflag=1
      return
    end if

    if (minval(abs(om_eigen)) == 0.0d0) then
      write(0,"(a)") "Error in inverse matrix operation. There exists a zero eigenvalue"
      errorflag=1
      return       
    end if

    allocate(om_inv(n,n), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in inverse matrix allocation"
      errorflag=1
      return
    end if
    om_inv=(0.0d0,0.0d0)

    do l=1,n
      do j=1,n
        do k=1,n
          if (abs(om_eigen(k))/maxval(abs(om_eigen)).gt.1.0d-3) then
            om_inv(l,j)=om_inv(l,j) + ovrlp(l,k)*dconjg(ovrlp(j,k))/om_eigen(k)
          end if
        end do
      end do
    end do

    do k=1,n
      do j=1,k-1
        if (dconjg(om_inv(k,j)).ne.om_inv(j,k)) then
          nonherm = nonherm + 1
        end if
      end do
    end do

    if (nonherm.ne.0) then
      write(0,"(a,i5,a,i5,a,i5)")"Error! Inverse matrix is Non-Hermitian. ",2*nonherm,&
              " unmatched off-diagonal elements found in a ", n,"x",n," matrix"
      errorflag=1
      return
    end if

    D = matmul(om_inv,C_k)

    deallocate(om_inv, stat = ierr)
    deallocate(om_eigen, stat = ierr)
    deallocate(WORK, stat = ierr)
    deallocate(RWORK, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of local arrays in matinv"
      errorflag=1
      return
    end if

    return

  end subroutine matinv2

!***********************************************************************************!

end module Ham    
