module clonecondense
  use alarrays
  use outputs
  use readpars

!*****************************************************************
!*   First attempt at translating the python file clonecombine.py into fortran to be called automatically
!*   and interact with bases much easier, as well as using pre-written functions.
!*
!*    Contains subroutines for:
!*    
!*     1) Recondensing the clones that are created through version 1 propagation
!*
!*
!*****************************************************************

contains 
  

  subroutine clone_condense(dt,timesteps, nbf, repeats, time_end,orgreps)
    implicit none

    type(basisfn), dimension(:), allocatable :: bs1, bs2
    character:: exdir, normpopfile !needs lengths
    character(len=12) clonedir
    character(len=4) repstr
    real, dimension(13) :: temparr
    real, allocatable,  dimension(:,:) :: normp, normweighting! combining all the norm files into this one matrix
    real, allocatable, dimension(:):: parent, child, time_clone, normWP, normWC, timestep_clone
    real, allocatable, dimension(:,:) :: shift
    integer(kind=4), intent(in):: nbf, timesteps, orgreps
    real(kind=8), intent(in):: dt, time_end
    real :: hold, crossterm1, crossterm2
    integer:: repeats, iostat, nlines, ierr, n, row1, row2, d, j, m

    ! Defining value of varibles
    ierr = 0
    nlines = 0
    crossterm1 = 0
    crossterm2 = 0
    ! Retrieving the working directory and finding the clone tag file
    ! ierr = getcwd(exdir)
    ! if (ierr.ne.0) stop 'getcwd:error'
    


    ! Opening the clone tag file and then reads in the different arrays
    clonedir = 'clonetag.out'
    open(11202,file=clonedir,status='old')
    
    DO
      READ(11202,*,iostat=ierr)
       IF (ierr/=0) EXIT
      nlines = nlines + 1
    END DO
    ! write(6,*) nlines
    allocate(parent(nlines))
    allocate(child(nlines))
    allocate(time_clone(nlines))
    allocate(timestep_clone(nlines))
    allocate(normWP(nlines))
    allocate(normWC(nlines))
    rewind(11202)
    do n=1, nlines
      read(11202,*) parent(n), child(n), time_clone(n), normWP(n), normWC(n)
    end do
    ! write(6,*) parent
    ! write(6,*) child
    ! write(6,*) time_clone
    
    timestep_clone = time_clone/time_end*timesteps
    write(6,*) timesteps
    ! write(6,*) timestep_clone
    ! write(6,*) normWP
    ! write(6,*) normWC
    
    ! Creates the norm weightings through information from the clone tag file.
    
    allocate(normweighting(timesteps+1,int(child(nlines))))
    allocate(normp(13,Timesteps+1))
    do n =1, int(child(nlines))
     do j =1, timesteps
      normweighting(j,n) = 1
     end do 
    end do
    do n =1, 13
      do j = 1, timesteps+1
       normp(n,j) = 0
      end do 
    end do

    ! write(6,*) size(normweighting), shape(normweighting)
    
    ! write(6,*) normweighting(1,:)
    do n=1, nlines
      ! write(6,*) 'here', n
      row1 = int(parent(n))
      ! write(6,*) row1
      row2 = int(child(n))
      ! write(6,*) row2
      hold = normweighting(int(timestep_clone(n)-1), row1)
      do m = int(timestep_clone(n)), timesteps+1
        ! write(6,*) 'here', n, m
        normweighting(m,row1) = normweighting(m,row1) * normWP(n)
        normweighting(m,row2) = hold * normWC(n)
      end do 
      do m =1, int(timestep_clone(n)-1)
        normweighting(m,row2) = 0
      end do
    end do 
    do n = 1, orgreps
      normweighting(2001,n) = normweighting(2000,n)
    end do
    
    ! write(6,*) normweighting(1,:)
    ! write(6,*) normweighting(1,:)
    ! write(6,*) normweighting(749,:)
    ! write(6,*) normweighting(750,:)
    ! write(6,*) normweighting(1499,:)
    ! write(6,*) normweighting(1500,:)
    write(6,*) normweighting(2000,:)
    write(6,*) normweighting(2001,:)

    ! Opening the different normpop files and collating them with their weighting attached. 

    do n=1, int(child(nlines))
      write(repstr,"(i4.4)") n
      normpopfile = 'normpop-'//trim(repstr)//'.out'
      write(6,*) 'normpop-'//trim(repstr)//'.out'
      open(11203, file = 'normpop-'//trim(repstr)//'.out')
      rewind(11203)
      read(11203,*)
      read(11203,*)
      read(11203,*)
      do j=1, timesteps+1
        read(11203,*) temparr(1), temparr(2), temparr(3), temparr(4), temparr(5), temparr(6), temparr(7), temparr(8), &
                      temparr(9), temparr(10), temparr(11), temparr(12), temparr(13)
        do m = 1, 13
          if (m.le.9) then
           normp(m,j) = normp(m,j) + temparr(m)
          end if
          if (m.ge.10) then
            normp(m,j) = normp(m,j) + temparr(m)*normweighting(j,n)
          end if 

        end do 
      end do 
      close(11203)
      ! write(6,*) normp
    end do

    ! Rescaling the weighted population matrix. 
    normp(1,:) = normp(1,:)/repeats
    normp(2,:) = normp(2,:)/repeats
    normp(3,:) = normp(3,:)/repeats
    normp(4,:) = normp(4,:)/repeats
    normp(5,:) = normp(5,:)/repeats
    normp(6,:) = normp(6,:)/repeats
    normp(7,:) = normp(7,:)/repeats
    normp(8,:) = normp(8,:)/repeats
    normp(9,:) = normp(9,:)/repeats
    normp(10,:) = normp(10,:)/orgreps
    normp(11,:) = normp(11,:)/orgreps
    normp(12,:) = normp(12,:)/orgreps
    normp(13,:) = normp(13,:)/orgreps
    write(6,*) normp(:,1)
    write(6,*) normp(:,100)
    write(6,*) normp(:,1000)
    write(6,*) normp(:,2000)
    write(6,*) normp(:,2001)

    ! Now to put that the population matrix into a collated output file.

    open(11204, file = 'normpop.out')
    write (11204,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
    write(11204,*)
    write(11204,*)
    do n = 1, timesteps+1
     write (11204,"(13(1x,es16.8e3))") normp(:,n)
    end do 
    close(11204)
  end subroutine

  subroutine cross_terms(filename1,filename2,normw1,normw2,nbf, crossterm1, crossterm2)
    type(basisfn), dimension (:), allocatable :: bs1, bs2
    real(kind=8), intent(inout) ::  crossterm1, crossterm2
    integer :: j, k, ierr
    complex(kind=8) :: sumamps1, sumamps2, ovrlp
    character(len=100), intent(in) :: filename1, filename2
    integer, intent(in) :: nbf
    real, intent(in) :: normw1, normw2



    sumamps1 = (0.0d0, 0.0d0)
    sumamps2 = (0.0d0, 0.0d0)
    crossterm1 = 0
    crossterm2 = 0

    ! Now to read in the basis sets themselves.
    call allocbs(bs1,nbf)
    ! call readcontbasis(bs1,trim(filename1),nbf)

    call allocbs(bs2,nbf)
    ! call readcontbasis(bs2,trim(filename2), nbf)

    !Calculating the overlap between the basis sets and summing over amplitudes

    do j=1,nbf
      do k=1,nbf
        ovrlp = ovrlp + ovrlpij(bs1(k)%z(:), bs2(j)%z(:))
        sumamps1 = sumamps1 + (dconjg(bs1(j)%a_pes(1))*bs2(k)%a_pes(1))
        sumamps2 = sumamps2 + (dconjg(bs1(j)%a_pes(2))*bs2(k)%a_pes(2))
      end do
    end do

    ! Calculating the cross terms

    crossterm1 = 2*real(sqrt((1-normw1)*normw2)*sumamps1*ovrlp)
    crossterm2 = 2*real(sqrt(normw1*(1-normw2))*sumamps2*ovrlp)


  end subroutine

    
end module