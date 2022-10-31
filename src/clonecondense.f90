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
    character:: exdir !needs lengths
    character(len=12) clonedir
    character(len=18) fn1, fn2
    character(len=4) repstr1, repstr2, repstr3, repstr4, repstr5
    real, dimension(13) :: temparr
    real, allocatable,  dimension(:,:) :: normp, normweighting ! combining all the norm files into this one matrix
    real, allocatable, dimension(:):: parent, child, time_clone, normWP, normWC, timestep_clone, rescale
    real, allocatable, dimension(:,:) :: shift, adjust1, adjust2
    integer(kind=4), intent(in):: timesteps, orgreps
    integer(kind=4), intent(inout):: nbf
    real(kind=8), intent(in):: dt, time_end
    real :: hold, crossterm1, crossterm2, timestep_parent, nw1,nw2, real_timestep
    integer:: repeats, iostat, nlines, ierr, n, row1, row2, d, j, m, place, jump, timestep, parent_rep

    ! Defining value of varibles
    ierr = 0
    nlines = 0
    crossterm1 = 0
    crossterm2 = 0
    jump = 1
    
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
    ! write(6,*) timesteps
    ! write(6,*) timestep_clone
    ! write(6,*) normWP
    ! write(6,*) normWC
    
    ! Allocates the shape of the matrices
    
    allocate(normweighting(timesteps+1,int(child(nlines))))
    allocate(normp(13,Timesteps+1))
    allocate(shift(repeats,3))
    allocate(adjust1(repeats,timesteps))
    allocate(adjust2(repeats,timesteps))
    allocate(rescale(timesteps))

    ! Initialising the matrices to their starting values before calculations

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

    do n = 1, 3
      do j = 1, repeats
        shift(n,j) = 0
      end do 
    end do 

    do n = 1, repeats
     do j = 1, timesteps
      adjust1(n,j) = 0 
      adjust2(n,j) = 0
     end do
    end do


    ! Creating the shift matrix 
    
    do n = 1, repeats
      shift(n,1) = n
      ! write(6,*) shift(n,1)
      if (any(n == child)) then 
        place = findloc(child, n, dim= 1)
        ! write(6,*) 'place for ', n, 'is', place
        ! write(6,*) timestep_clone(place), parent(place)
        shift(n,2) = timestep_clone(place)
        shift(n,3) = parent(place)
        
      end if 
    end do 
    do n = 1, repeats
      write(6,*) shift(n,:)
    end do 



    ! Creates the norm weightings through information from the clone tag file.
    
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
    ! write(6,*) normweighting(2000,:)
    ! write(6,*) normweighting(2001,:)


     ! creating the 2 adjust matrices (comment out if no crossterms)

    do n = 1, repeats 
      adjust1(n,:) = 0
      adjust2(n,:) = 0 
    end do 
    do n = repeats, orgreps + 1, -1
      m=1 
      real_timestep = 0
      write(6,*) 'n, ', n
      
      
      do while (real_timestep.lt.timesteps)
        ! write(6,*) 'm, ', m 
        real_timestep = shift(n,2) + (m*jump) 
        parent_rep = shift(n,3)
        timestep_parent = real_timestep - shift(parent_rep,2)
        write(6,*) 'real_timestep, ', real_timestep
        write(repstr1,"(i4.4)") n
        write(repstr2,"(i4.4)") int(shift(n,3))
        write(repstr3,"(i4.4)") int(m*jump)
        write(repstr4,"(i4.4)") int(timestep_parent)
        
        fn1 = 'outbscon-'//repstr1//'-'//repstr3 
        fn2 = 'outbscon-'//repstr2//'-'//repstr4
        write(6,*) fn1, fn2
        nw1 = normweighting(int(real_timestep),n)
        nw2 = normweighting(int(real_timestep),int(shift(n,3)))
        write(6,*) 'normweighting, ', nw1,nw2
        call cross_terms(fn1,fn2,nw1,nw2, nbf, crossterm1,crossterm2)
        ! write(6,*) 'crossterms, ', crossterm1, crossterm2
        write(6,*) '****************'
        do j = int(timestep_clone(n)+(m*jump)), timesteps
          adjust1(n,j) = crossterm1
          adjust2(n,j) = crossterm2
        end do
        write(6,*) 'crossterm1 is, ', crossterm1, 'and so adjust1(n,j) = ', adjust1(n,int(timestep_clone(n)+(m*jump)))
        if ((timestep_clone(n)+m*jump).eq.timesteps) then
          call cross_terms(fn1,fn2,nw1,nw2, nbf, crossterm1,crossterm2)
          adjust1(n,timesteps) = crossterm1
          adjust2(n,timesteps) = crossterm2
        else 
          adjust1(n,timesteps) = adjust1(n,timesteps-1)
          adjust2(n,timesteps) = adjust2(n,timesteps-1)
        end if 
        m = m+1
      end do  
    end do
    write(6,*) 'adjust1 at 0, ', adjust1(:,0)
    write(6,*) 'adjust1 at 751, ', adjust1(:,751)
    write(6,*) 'adjust2 at 0, ', adjust2(:,0) 
    write(6,*) 'adjust2 at 751, ', adjust2(:,751)

    ! Opening the different normpop files and collating them with their weighting attached. 

    do n=1, int(child(nlines))
      write(repstr5,"(i4.4)") n
      write(6,*) 'normpop-'//trim(repstr5)//'.out'
      open(11203, file = 'normpop-'//trim(repstr5)//'.out')
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
    ! write(6,*) normp(:,1)
    ! write(6,*) normp(:,100)
    ! write(6,*) normp(:,1000)
    ! write(6,*) normp(:,2000)
    ! write(6,*) normp(:,2001)

    ! Time to use the adjust matrices to, well, adjust. (comment this out if no crossterms)

    do n=1, repeats
     do m = 1, timesteps
       normp(10,m) = normp(10,m) + adjust1(n,m)
       normp(11,m) = normp(11,m) + adjust2(n,m)
      end do 
    end do 
    
    do m = 1, timesteps
      rescale(m) = normp(10,m)+normp(11,m)
      normp(10,m) = normp(10,m)/rescale(m)
      normp(11,m) = normp(11,m)/rescale(m)
      normp(12,m) = normp(10,m)+normp(11,m)
      normp(13,m) = normp(10,m)-normp(11,m)
    end do 




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
    real, intent(inout) ::  crossterm1, crossterm2
    integer :: j, k, ierr
    complex(kind=8), dimension(nbf,nbf) ::ovrlp12, ovrlp21
    complex(kind=8) :: sumamps1, sumamps2, ovrlp1, ovrlp2, cpop11, cpop12, cpop21,cpop22
    character(len=18), intent(in) :: filename1, filename2
    integer, intent(inout) :: nbf
    real, intent(in) :: normw1, normw2



    sumamps1 = (0.0d0, 0.0d0)
    sumamps2 = (0.0d0, 0.0d0)
    crossterm1 = 0
    crossterm2 = 0
    ovrlp1 = (0.0d0,0.0d0)
    ovrlp2 = (0.d0,0.0d0)
  

    ! Now to read in the basis sets themselves.
    call allocbs(bs1,nbf)
    call readcontbasis(bs1,trim(filename1),nbf)
    

    call allocbs(bs2,nbf)
    call readcontbasis(bs2,trim(filename2), nbf)
    
  
    !Calculating the overlap between the basis sets and summing over amplitudes
    ovrlp12 = ovrlpmat12(bs1,bs2)
    ovrlp21 = ovrlpmat21(bs1,bs2)
    
    ! FIRST ATTEMPT, POSSIBLE GARBAGE
    ! do j=1,nbf
    !   do k=j,nbf
    !     ! write(6,*) ovrlp
    !     sumamps1 = sumamps1 + (dconjg(bs1(j)%a_pes(1))*bs2(k)%a_pes(1))
    !     ! write(6,*) sumamps1
    !     sumamps2 = sumamps2 + (dconjg(bs1(j)%a_pes(2))*bs2(k)%a_pes(2))
    !     ! write(6,*) sumamps2
    !   end do
    !   do k=j,nbf
    !     ovrlp1 = ovrlp1 + ovrlpij(bs1(k)%z, bs2(j)%z) * sumamps1
    !     ovrlp2 = ovrlp2 + ovrlpij(bs2(k)%z, bs1(j)%z) * sumamps2
    !   end do 
    ! end do
    
    ! write(6,*) ovrlp1, ovrlp2, sumamps1,sumamps2
    ! write(6,*) ovrlp1*sumamps1
    ! write(6,*) ovrlp2*sumamps2
    ! write(6,*) normw1, normw2
    ! crossterm1 = 2*real(sqrt((1-normw1)*normw2)*ovrlp1)
    ! crossterm2 = 2*real(sqrt(normw1*(1-normw2))*ovrlp2)
    ! Calculating the cross terms

    do k=1,nbf
      do j=1,nbf
        cpop11 = cpop11 + dconjg(bs1(j)%D_big) * ovrlp12(j,k) * bs2(k)%D_big &
                * dconjg(bs1(j)%d_pes(1)) * bs2(k)%d_pes(1) &
                * cdexp(i*(bs1(k)%s_pes(1)-bs2(j)%s_pes(1)))
      end do
    end do
    do k=1,nbf
      do j=1,nbf
        cpop12 = cpop12 + dconjg(bs2(j)%D_big) * ovrlp21(j,k) * bs1(k)%D_big &
                * dconjg(bs2(j)%d_pes(1)) * bs1(k)%d_pes(1) &
                * cdexp(i*(bs2(k)%s_pes(1)-bs1(j)%s_pes(1)))
      end do
    end do

    do k=1,nbf
      do j=1,nbf
        cpop21 = cpop21 + dconjg(bs1(j)%D_big) * ovrlp12(j,k) * bs2(k)%D_big &
                * dconjg(bs1(j)%d_pes(2)) * bs2(k)%d_pes(2) &
                * cdexp(i*(bs1(k)%s_pes(2)-bs2(j)%s_pes(2)))
      end do
    end do

    do k=1,nbf
      do j=1,nbf
        cpop22 = cpop22 + dconjg(bs2(j)%D_big) * ovrlp21(j,k) * bs1(k)%D_big &
                * dconjg(bs2(j)%d_pes(2)) * bs1(k)%d_pes(2) &
                * cdexp(i*(bs2(k)%s_pes(2)-bs1(j)%s_pes(2)))
      end do
    end do

    crossterm1 = sqrt(normw1*normw2)*(cpop11+cpop12) 
    crossterm2 = sqrt(normw1*normw2)*(cpop21+cpop22)
    !
    write(6,*) 'crossterms, ', crossterm1, crossterm2


  end subroutine

  function ovrlpmat12(bs1,bs2)

    implicit none
    type(basisfn), dimension (:), intent(in)::bs1, bs2
    complex(kind=8), dimension(size(bs1),size(bs1))::ovrlpmat12
    integer::k,j,n

    if (errorflag .ne. 0) return

    ovrlpmat12 = (0.0d0, 0.0d0)

    n = size(bs1)

    do k=1,n
      do j=k,n
        if (j==k) then
          ovrlpmat12(j,k)=(1.0d0,0.0d0)
        else
          ovrlpmat12(j,k)=ovrlpij(bs1(j)%z,bs2(k)%z)
          ovrlpmat12(k,j)=dconjg(ovrlpmat12(j,k))
        end if
      end do
    end do
    
    return

  end function ovrlpmat12

  function ovrlpmat21(bs1,bs2)

    implicit none
    type(basisfn), dimension (:), intent(in)::bs1, bs2
    complex(kind=8), dimension(size(bs1),size(bs1))::ovrlpmat21
    integer::k,j,n

    if (errorflag .ne. 0) return

    ovrlpmat21 = (0.0d0, 0.0d0)

    n = size(bs1)

    do k=1,n
      do j=k,n
        if (j==k) then
          ovrlpmat21(j,k)=(1.0d0,0.0d0)
        else
          ovrlpmat21(j,k)=ovrlpij(bs2(j)%z,bs1(k)%z)
          ovrlpmat21(k,j)=dconjg(ovrlpmat21(j,k))
        end if
      end do
    end do
    
    return

  end function ovrlpmat21
    
end module