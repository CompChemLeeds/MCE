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

!***********************************************************************************!
!            Hamiltonian Subroutines and Functions
!***********************************************************************************!

  function Hehr(bf,t,reps)   !   Level 1 Function

    implicit none
    type (basisfn), intent(inout)::bf
    real(kind=8), intent (in) :: t
    integer, intent(in) :: reps
    
    complex(kind=8)::Hehr, asum
    complex(kind=8),dimension(:,:), allocatable :: z
    complex(kind=8),dimension(:,:,:), allocatable::H
    integer::m,r,s,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(H(1,npes,npes), stat = ierr)
    if (ierr==0) allocate (z(1,ndim), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of H matrix in Hehr"
      errorflag=1
      return
    end if
    
    do m=1,ndim
      z(1,m)=bf%z(m)
    end do

    Hehr = (0.0d0, 0.0d0)
    asum = (0.0d0, 0.0d0)

    call Hijdiag(H,z,t,reps) 

    do r=1,npes
      do s=1,npes
        Hehr = Hehr + H(1,r,s)*dconjg(bf%a_pes(r))*bf%a_pes(s)
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

!***********************************************************************************!

end module Ham    
