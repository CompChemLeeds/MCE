MODULE alarrays

  use globvars

!***********************************************************************************!
!*
!*         Allocation/Deallocation Module
!*           
!*   Contains subroutines for:
!*
!*      1) Allocating & initialising the basis set
!*      2) Allocating & initialising a single basis function
!*      3) Deallocating a single basis function
!*      4) Deallocating the basis set
!*      5) Allocating and initialising the multi-configurational nbf x nbf 
!*           Hamiltonian matrix
!*      6) Allocating and initialising the single configurational npes x npes 
!*           Hamiltonian matrices
!*      7) Deallocating the multi-configurational Hamiltonian matrix
!*      8) Deallocating the single-configurational Hamiltonian matrices
!*      
!***********************************************************************************!

contains

!***********************************************************************************!
!               Array Allocation/Deallocation Subroutines
!***********************************************************************************!

  subroutine allocbs(bs, nbf)
  
    ! Allocates the top level of the variable of the basisfn defined type

    implicit none
    
    type(basisfn),dimension(:),allocatable, intent(inout)::bs
    integer, intent (in) :: nbf
    
    integer:: j, ierr

    if (errorflag .ne. 0) return

    ierr=0
    if (allocated(bs).eqv..false.) then
      allocate (bs(nbf), stat=ierr) !allocate the top level of the defined type array 
      if (ierr/=0) then
        write(0,"(a,i0)") "Error in basis set allocation. ierr had value ", ierr
        errorflag=1
        return
      end if
    end if

    do j=1,size(bs)
      call allocbf(bs(j))! this allocates each bf within bs
    end do

    return

  end subroutine allocbs

!------------------------------------------------------------------------------------

  subroutine allocbf(bf)
  
    ! Allocates the elements of the variable of the basisfn defined type

    implicit none
    
    type(basisfn),intent(inout)::bf
    
    integer :: ierr

    if (errorflag .ne. 0) return

    ierr=0

    ! Each sub-array must be allocated separately. 
    ! Allocation continues only if no errors occur in previous allocations

    allocate (bf%z(ndim), stat=ierr)                     
    if (ierr==0) allocate (bf%d_pes(npes), stat=ierr)    
    if (ierr==0) allocate (bf%s_pes(npes), stat=ierr)    
    if (ierr==0) allocate (bf%a_pes(npes), stat=ierr)  
    if (ierr==0) allocate (bf%carray(5), stat=ierr)  
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in basis function allocation. ierr had value ", ierr
      errorflag=1
      return
    end if


    bf%z    (1:ndim) = (0.0d0,0.0d0)
    bf%d_pes(1:npes) = (0.0d0,0.0d0)
    bf%a_pes(1:npes) = (0.0d0,0.0d0)
    bf%s_pes(1:npes) =  0.0d0
    bf%D_big         = (0.0d0,0.0d0)
    bf%carray(1:5)  =  0.0d0

    return

  end subroutine allocbf

!------------------------------------------------------------------------------------

  subroutine deallocbf(bf)
  
    ! Deallocates the elements of the variable of the basisfn defined type

    implicit none
    
    type(basisfn),intent(inout)::bf
    
    integer :: ierr

    if (errorflag .ne. 0) return

    ierr=0

    ! Each sub-array must be deallocated separately.
    ! Deallocation continues only if no errors occur in previous allocations

    deallocate (bf%z, stat=ierr)
    if (ierr==0) deallocate (bf%d_pes, stat=ierr)
    if (ierr==0) deallocate (bf%s_pes, stat=ierr)
    if (ierr==0) deallocate (bf%a_pes, stat=ierr)
    if (ierr==0) deallocate (bf%carray, stat=ierr)
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in basis function deallocation. ierr had value ", ierr
      errorflag=1
      return
    end if


    return

  end subroutine deallocbf

!------------------------------------------------------------------------------------

  subroutine deallocbs(bs)
  
    ! Deallocates the top level of the variable of the basisfn defined type

    implicit none
    
    type(basisfn),dimension(:), allocatable, intent(inout)::bs
    
    integer:: j, ierr

    if (errorflag .ne. 0) return

    ierr=0

    do j=1,size(bs)
      call deallocbf(bs(j))! this deallocates each bf within bs
    end do

    ! Deallocation is carried out bottom up. 
    ! Top level deallocation only occurs if lower level deallocation was successful

    if (errorflag==0) deallocate (bs, stat=ierr)
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in basis set deallocation. ierr had value ", ierr
      errorflag=1
      return
    end if

    return

  end subroutine deallocbs

!------------------------------------------------------------------------------------

  subroutine allocham(H, nbf)
  
    ! Allocates the top level of the variable of the hamiltonian defined type

    implicit none
    
    type(hamiltonian),dimension(:,:),allocatable,intent (inout)::H
    integer, intent (in) :: nbf
    
    integer::k, j, ierr

    if (errorflag .ne. 0) return

    ierr=0

    ! Allocate the top level of the array before calling the subroutine to allocate
    ! the different sub-arrays

    allocate (H(nbf, nbf), stat=ierr)
    do k=1,nbf
      do j=1,nbf
        call allochamjk(H(j,k))
      end do
    end do
    if (ierr/=0) then
      write(0,"(2a,i0)") "Error in hamiltonian matrix allocation.", &
                          " ierr had value ", ierr
      errorflag=1
      return
    end if

  end subroutine allocham

!------------------------------------------------------------------------------------

  subroutine allochamjk(Hsing)
  
    ! Allocates the elements of the variable of the hamiltonian defined type

    implicit none
    
    type(hamiltonian),intent(inout) :: Hsing
    
    integer :: ierr
    
    if (errorflag .ne. 0) return

    ierr = 0
    
    ! Each sub-array has to be allocated separately

    allocate (Hsing%Hjk(npes, npes), stat=ierr)
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in hamiltonian element allocation.",&
                          " ierr had value ", ierr
      errorflag=1
      return
    end if
    Hsing%Hjk = (0.0d0,0.0d0)

    return

  end subroutine allochamjk

!------------------------------------------------------------------------------------

  subroutine deallocham(H)
  
    ! Deallocates the top level of the variable of the hamiltonian defined type

    implicit none
    
    type(hamiltonian),dimension(:,:),allocatable,intent (inout)::H
    
    integer::k, j, ierr

    if (errorflag .ne. 0) return

    ierr=0
    
    ! Deallocate the top level of the array before calling the subroutine to 
    ! deallocate the different sub-arrays

    do k=1,size(H,2)
      do j=1,size(H,1)
        call deallochamjk(H(j,k))
      end do
    end do
    deallocate (H, stat=ierr)
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in hamiltonian matrix deallocation.",&
                          " ierr had value ", ierr
      errorflag=1
      return
    end if

  end subroutine deallocham

!------------------------------------------------------------------------------------

  subroutine deallochamjk(Hsing)
  
    ! Deallocates the elements of the variable of the hamiltonian defined type

    implicit none
    
    type(hamiltonian),intent(inout) :: Hsing
    
    integer :: ierr
    
    if (errorflag .ne. 0) return

    ierr = 0
    
    ! Each sub-array is deallocated separately

    deallocate (Hsing%Hjk, stat=ierr)
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in hamiltonian element deallocation.",&
                          " ierr had value ", ierr
      errorflag=1
      return
    end if

    return

  end subroutine deallochamjk

!***********************************************************************************!

END MODULE alarrays
