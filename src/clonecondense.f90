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
  

  subroutine clone_condense(dt,timesteps, nbf, repeats, time_end,orgreps, clonefreq)
    implicit none

    type(basisfn), dimension(:), allocatable :: bs1, bs2
    real, dimension(:,:,:), allocatable :: populations, ctarray, normpfs, condreps
    character:: exdir !needs lengths
    character(len=12) clonedir
    character(len=18) fn
    character(len=4) repstr1, repstr2, repstr3, repstr4, repstr5
    real, dimension(13) :: temparr
    real, allocatable,  dimension(:,:) :: normp_final ! combining all the norm files into this one matrix
    real, allocatable, dimension(:):: parent, child, time_clone, initrep, timestep_clone, clone_event,  normWP, normWC 
    real, allocatable, dimension(:,:) :: shift, normweighting
    integer, dimension(:,:), allocatable :: tree
    integer,dimension(:), allocatable :: normwplace
    integer(kind=4), intent(in):: timesteps, orgreps, clonefreq
    integer(kind=4), intent(inout):: nbf
    real(kind=8)::crossterm1, crossterm2
    real(kind=8), intent(in):: dt, time_end
    integer, dimension(:), allocatable :: index0, farray, clonestep, tree_hold, index_clone 
    real :: hold,  timestep_parent, nw1,nw2, real_timestep, ct1, ct2, num_events, rescale
    integer:: repeats, iostat, nlines, ierr, n, row1, row2, d, j, m, place, jump, timestep, parent_rep, k, l, npunit, npunit2
    integer:: h, f, i, e, ix, normflg, i1,i2,i3,i4



    ! Defining value of varibles
    ierr = 0
    nlines = 0
    crossterm1 = 0
    crossterm2 = 0
    jump = 1
    
    
    ! Retrieving the working directory and finding the clone tag file
    ! ierr = getcwd(exdir)
    ! if (ierr.ne.0) stop 'getcwd:error'
    
    ! FLAGS FOR METHODS (0 for mine, 1 for Oliver's)
    normflg = 1.0d0

    ! Opening the clone tag file and then reads in the different arrays
    clonedir = 'clonetag.out'
    open(11202,file=clonedir,status='old')

    DO
      READ(11202,*,iostat=ierr)
       IF (ierr/=0) EXIT
      nlines = nlines + 1
    END DO
    write(6,*) nlines

    ! write(6,*) nlines
    allocate(parent(nlines))
    allocate(child(nlines))
    allocate(time_clone(nlines))
    allocate(timestep_clone(nlines))
    allocate(initrep(nlines))
    allocate(normWC(nlines))
    allocate(normWP(nlines))
    allocate(populations(nlines,timesteps,2))
    allocate(ctarray(repeats,timesteps,2))
    allocate(normpfs(13,timesteps,repeats))
    allocate(condreps(orgreps,13,timesteps))
    
    allocate(normwplace(repeats))

    if (mod(timesteps,clonefreq).eq.0) then
      num_events = timesteps/clonefreq - 1 
    else 
      num_events = timesteps/clonefreq 
    end if 
    write(6,*) 'number of cloning events is, ', num_events
    allocate(clone_event(num_events+1))
    allocate(normweighting(repeats,num_events+1))
    write(6,*) 'size of nw = ,', size(normweighting)
    
    do n =1, size(clone_event)-1
      clone_event(n) = n*clonefreq
    end do
    clone_event(num_events+1) = timesteps 
    write(6,*) clone_event
    
    
   
    
    

    rewind(11202)
    do n=1, nlines
      read(11202,*) parent(n), child(n), time_clone(n), initrep(n), normWP(n), normWC(n)
    end do
    ! write(6,*) parent
    ! write(6,*) child
    ! write(6,*) time_clone
    ! write(6,*) initrep
    ! write(6,*) normWP
    ! write(6,*) normWC
    close(11202)
    write(6,*) '*********clonetag read********'
    timestep_clone = time_clone/time_end*timesteps
    ! write(6,*) timesteps
    ! write(6,*) timestep_clone
   
    
    ! Allocates the shape of the matrices
    
  
    allocate(normp_final(13,timesteps+1))
    allocate(shift(repeats,3))
    allocate(tree(orgreps, 2**num_events))
    allocate(tree_hold(2**num_events))
    allocate(clonestep(2**num_events))
    allocate(farray(2**num_events))
    

    ! Initialising the matrices to their starting values before calculations

    tree = 0.d0
    shift = 0.d0
    normp_final = 0.0d0
    farray = 0.d00
    condreps = 0.d00
    ctarray = 0.0d0
    normweighting = 1.0d0
    normwplace = 2.0d0
    
    ! Creating the shift matrix 
   
    do n = 1, repeats
      shift(n,1) = n
      ! write(6,*) shift(n,1)
      if (any(n == child)) then 
        index0 = minloc((child-n)**2)
        ! write(6,*) 'back up place for ', n, 'is', index0(1)
        ! write(6,*) timestep_clone(index0(1)), parent(index0(1))
        shift(n,2) = timestep_clone(index0(1))
        shift(n,3) = parent(index0(1))
      end if 
    end do 
    do n = 1, repeats
      write(6,*) shift(n,:)
    end do 
    write(6,*) '******** shift matrix made ********'
    
    ! Creating the tree array

    do n=1, orgreps
      tree(n,1) = n
    end do 
    do n=1, nlines
      farray = 0.d0
      ! write(6,*) child(n), initrep(n) 
      do k = 1, size(farray)
        farray(k) = tree(initrep(n), k)
      end do
      m = n+ orgreps
      ! write(6,*) 'fake array = ', farray
      index0 = minloc(farray)
      ! write(6,*) 'here index = ', initrep(n), index0(1), ' and input is = ', child (n)
      
      tree(initrep(n), index0(1)) = child(n)
    end do 


    ! do n=1,orgreps
    !  write(6,*) tree(n,:)
    ! end do

    write(6,*) '***** Tree matrix made ******'
    
    
    do f=1,orgreps
      tree_hold = 0
      tree_hold(1) = tree(f,1)
      e = 1 
      k=1
      do m=1,2**num_events
        clonestep(m) = shift(tree(f,m),2)
      end do 
      ! write(6,*) 'cloning events ', clonestep
      do n=1,size(clone_event)-1
        ! write(6,*) 'clonestep = ', clone_event(n)
        
        index_clone = pack([(ix,ix=1,size(clonestep))],clonestep==clone_event(n))
        ! write(6,*) 'indices for the clone  ', index_clone
        
        index0 = minloc(tree_hold)
        ! write(6,*) 'indices for temp array ', index0
        do m = 1, size(index_clone)
          tree_hold(index0(1)+(m-1)) = tree(f,index_clone(m))
        end do
        e = e + k
        k=k*2
        
      end do 
      tree(f,:) = tree_hold(:)

    end do 
    do n=1,orgreps
     write(6,*) tree(n,:)
    end do

    write(6,*) '***** Tree matrix ordered ******'
    
    ! Creates the norm weightings through information from the clone tag file.
    
  

   

    

    normweighting = 1.0d0
    do n = 1, nlines
      row1 = int(parent(n))
      write(6,*) 'row1 is ', row1
      row2 = int(child(n))
      write(6,*) 'row 2 is, ', row2
      normweighting(row2, normwplace(row1)) = normweighting(row1,normwplace(row1)-1) * normWC(n)
      normweighting(row1, normwplace(row1)) = normweighting(row1,normwplace(row1)-1) * normWP(n)
      normwplace(row1) = normwplace(row1) + 1
      normwplace(row2) = normwplace(row1)
    end do 

    do l=1, repeats
      write(6,*) normweighting(l,:)
    end do
    
    
    
    
   
    

  
   
    
    write(6,*) '***** Normweighting matrix made *****' 

    do l=1, repeats
      write(repstr1,"(i4.4)") l
      write(6,*) repstr1
      open(11203, file = 'normpop-'//trim(repstr1)//'.out')
      read(11203,*)
      read(11203,*)
      read(11203,*)
      do j = 1, timesteps
        read(11203,*) temparr(1), temparr(2), temparr(3), temparr(4), temparr(5), temparr(6), temparr(7), temparr(8), &
                      temparr(9), temparr(10), temparr(11), temparr(12), temparr(13)
        do m = 1, 13
         normpfs(m,j,l) =  normpfs(m,j,l) + temparr(m)
        end do 
        ! write(6,*) 'rep ', l, 'timestep ', j, 'populations ', normpfs(10,j,l), normpfs(11,j,l)
      end do 
      close(11203)
      write(6,*) 'Read in file, ', l  
    end do 


    n=1
    do k=1, size(clone_event)-1
      n = n*2
      do j=clone_event(k), clone_event(k+1)-1
        do i=1,orgreps
          ! write(6,*) '********i = ', i, '*********'
          do h=1, n    
            do f = h+1, n 
              i1 = tree(i,h)
              i2 = int(j-shift(tree(i,h),2))
              i3 = tree(i,f)
              i4 = int(j-shift(tree(i,f),2))
              ! write(6,*) tree(i,h), int(j-shift(tree(i,h),2)), tree(i,f), int(j-shift(tree(i,f),2))
              ! write(6,*) fn1, ' ', fn2
              if (normflg==1) then 
                if (j==clone_event(k)) then 
                  nw1 = normweighting(tree(i,h),k+1)
                  nw2 = normweighting(tree(i,f),k+1)
                else 
                  nw1 = normweighting(tree(i,h),k+1)
                  nw2 = normweighting(tree(i,f),k+1)
                end if
              else 
                nw1 = 1.0d0
                nw2 = 1.0d0
              end if 
              call cross_terms(i1,i2,i3,i4,nbf, crossterm1,crossterm2,nw1,nw2)
              ! write(6,*) crossterm1, crossterm2
              ! write(6,*) '****************'
              ctarray(i,j,1) = ctarray(i,j,1) + crossterm1
              ctarray(i,j,2) = ctarray(i,j,2) + crossterm2
              ! write(6,*) ctarray(i,j,1), ctarray(i,j,2)
            end do
          end do
          ! write(6,*) 'crossterms, ', i, j, ctarray(i,j,1) , ctarray(i,j,2)
        end do
        
      end do
    end do
    
    

    !!$omp parallel
    !!$omp do  &
    !!$omp private(i,j,m,n,rescale) &
    !!$omp shared(condreps,normweighting,normpfs,tree,normp_final)
    ! do i =1, orgreps ! 1
    !   do j=1, clone_event(1) !2
    !     do m =1, 13 !6 
    !       condreps(i,m,j) = normpfs(m,j,i)
    !     end do !6
    !   end do !2
    !   e=2
    !   do n=1, size(clone_event)-1 !3
    !     ! write(6,*) clone_event(n)+1, clone_event(n+1)
    !     do j=clone_event(n)+1, clone_event(n+1) !4
    !       do f=1,e !5
    !         condreps(i,1,j) = condreps(i,1,j) + normweighting(tree(i,f),n+1)*normpfs(1,j,tree(i,f))
    !         condreps(i,2,j) = condreps(i,2,j) + normweighting(tree(i,f),n+1)*normpfs(2,j,tree(i,f))
    !         condreps(i,3,j) = condreps(i,3,j) + normweighting(tree(i,f),n+1)*normpfs(3,j,tree(i,f))
    !         condreps(i,4,j) = condreps(i,4,j) + normweighting(tree(i,f),n+1)*normpfs(4,j,tree(i,f))
    !         condreps(i,5,j) = condreps(i,5,j) + normweighting(tree(i,f),n+1)*normpfs(5,j,tree(i,f))
    !         condreps(i,6,j) = condreps(i,6,j) + normweighting(tree(i,f),n+1)*normpfs(6,j,tree(i,f))
    !         condreps(i,7,j) = condreps(i,7,j) + normweighting(tree(i,f),n+1)*normpfs(7,j,tree(i,f))
    !         condreps(i,8,j) = condreps(i,8,j) + normweighting(tree(i,f),n+1)*normpfs(8,j,tree(i,f))
    !         condreps(i,9,j) = condreps(i,9,j) + normweighting(tree(i,f),n+1)*normpfs(9,j,tree(i,f))
    !         condreps(i,10,j) = condreps(i,10,j) + normweighting(tree(i,f),n+1)*normpfs(10,j,tree(i,f))
    !         condreps(i,11,j) = condreps(i,11,j) + normweighting(tree(i,f),n+1)*normpfs(11,j,tree(i,f))
    !       end do ! 5
    !       write(6,*) 'pops 1&2 with ct 1&2', condreps(i,10,j), condreps(i,11,j), ctarray(i,j,1), ctarray(i,j,2)
    !       condreps(i,10,j) = condreps(i,10,j) + ctarray(i,j,1)
    !       condreps(i,11,j) = condreps(i,11,j) + ctarray(i,j,2) 
    !       write(6,*) 't =, ', j, 'tot pops with crossterms ', condreps(i,10,j), condreps(i,11,j)
    !       rescale = condreps(i,10,j) + condreps(i,11,j)
    !       condreps(i,10,j) = condreps(i,10,j)/rescale
    !       condreps(i,11,j) = condreps(i,11,j)/rescale
    !       condreps(i,12,j) = condreps(i,10,j) + condreps(i,11,j)
    !       condreps(i,13,j) = condreps(i,10,j) - condreps(i,11,j)
    !       write(6,*) 'tot pops with crossterms rescaled', condreps(i,10,j), condreps(i,11,j)
    !       write(6,*) 'orgrep =, ', i, 't = ', j, 'p1+p2 ', condreps(i,12,j), 'p1-p2 ', condreps(i,13,j)
    !       write(6,*) '********************'
    !     end do ! 4
    !     e =e*2
    !   end do !3
    ! end do !1
    !!$omp end do

    !!$omp do collapse(2)


    do i =1, orgreps ! 1
      do j=1, clone_event(1) !2
        do m =1, 13 !6 
          condreps(i,m,j) = normpfs(m,j,i)
        end do !6
      end do !2
      e=2
      do n=1, size(clone_event)-1 !3
        ! write(6,*) clone_event(n)+1, clone_event(n+1)
        do j=clone_event(n)+1, clone_event(n+1) !4
          do f=1,e !5
            condreps(i,1,j) = condreps(i,1,j) + normweighting(tree(i,f),n+1)*normpfs(1,j,tree(i,f))
            condreps(i,2,j) = condreps(i,2,j) + normweighting(tree(i,f),n+1)*normpfs(2,j,tree(i,f))
            condreps(i,3,j) = condreps(i,3,j) + normweighting(tree(i,f),n+1)*normpfs(3,j,tree(i,f))
            condreps(i,4,j) = condreps(i,4,j) + normweighting(tree(i,f),n+1)*normpfs(4,j,tree(i,f))
            condreps(i,5,j) = condreps(i,5,j) + normweighting(tree(i,f),n+1)*normpfs(5,j,tree(i,f))
            condreps(i,6,j) = condreps(i,6,j) + normweighting(tree(i,f),n+1)*normpfs(6,j,tree(i,f))
            condreps(i,7,j) = condreps(i,7,j) + normweighting(tree(i,f),n+1)*normpfs(7,j,tree(i,f))
            condreps(i,8,j) = condreps(i,8,j) + normweighting(tree(i,f),n+1)*normpfs(8,j,tree(i,f))
            condreps(i,9,j) = condreps(i,9,j) + normweighting(tree(i,f),n+1)*normpfs(9,j,tree(i,f))
            condreps(i,10,j) = condreps(i,10,j) + normweighting(tree(i,f),n+1)*normpfs(10,j,tree(i,f))
            condreps(i,11,j) = condreps(i,11,j) + normweighting(tree(i,f),n+1)*normpfs(11,j,tree(i,f))
          end do ! 5
          !write(6,*) 'pops 1&2 with ct 1&2', condreps(i,10,j), condreps(i,11,j), ctarray(i,j,1), ctarray(i,j,2)
          condreps(i,10,j) = condreps(i,10,j) + ctarray(i,j,1)
          condreps(i,11,j) = condreps(i,11,j) + ctarray(i,j,2) 
          !write(6,*) 't =, ', j, 'tot pops with crossterms ', condreps(i,10,j), condreps(i,11,j)
          rescale = condreps(i,10,j) + condreps(i,11,j)
          condreps(i,10,j) = condreps(i,10,j)/rescale
          condreps(i,11,j) = condreps(i,11,j)/rescale
          condreps(i,12,j) = condreps(i,10,j) + condreps(i,11,j)
          condreps(i,13,j) = condreps(i,10,j) - condreps(i,11,j)
          ! write(6,*) 'tot pops with crossterms rescaled', condreps(i,10,j), condreps(i,11,j)
          ! write(6,*) 'orgrep =, ', i, 't = ', j, 'p1+p2 ', condreps(i,12,j), 'p1-p2 ', condreps(i,13,j)
          ! write(6,*) '********************'
        end do ! 4
        e =e*2
      end do !3
    end do !1
    do i =1, orgreps
      do j=1, timesteps
        do m=1,13
          normp_final(m,j) = normp_final(m,j) + condreps(i,m,j)/orgreps
        end do 
      end do
    end do
    !!$omp end do 

  !!$omp end parallel 
      

    ! Now to put that the population matrix into a collated output file.

    open(11204, file = 'normpop.out')
    write (11204,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
    write(11204,*)
    write(11204,*)
    do n = 1, timesteps
     write (11204,"(13(1x,es16.8e3))") normp_final(:,n)
    end do 
    close(11204)

    do i =1, orgreps
      write(repstr1,"(i4.4)") i 
      ! npunit = 1020 + i 
      open(11203, file = 'condnorm-'//trim(repstr1)//'.out')
      write (11203,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
      write(11203,*)
      write(11203,*)
      do n = 1, timesteps
        write (11203,"(13(1x,es16.8e3))") condreps(i,:,n)
      end do 
      close(11203)
    end do 


    write(6,*) 'CONDENSING FINISHED'
    
  end subroutine







  subroutine cross_terms(i1,i2,i3,i4,nbf, crossterm1, crossterm2,nw1,nw2)
    type(basisfn), dimension (:), allocatable :: bs1, bs2
    real(kind=8), intent(inout) ::  crossterm1, crossterm2
    real(kind=4), intent(in):: nw1, nw2
    integer, intent(in):: i1,i2,i3,i4
    integer :: j, k, ierr
    complex(kind=8), dimension(nbf,nbf) ::ovrlp12, ovrlp21
    complex(kind=8) :: sumamps1, sumamps2, ovrlp1, ovrlp2, cpop11, cpop12, cpop21,cpop22
    integer, intent(inout) :: nbf
    



    sumamps1 = (0.0d0, 0.0d0)
    sumamps2 = (0.0d0, 0.0d0)
    crossterm1 = 0
    crossterm2 = 0
    ovrlp1 = (0.0d0,0.0d0)
    ovrlp2 = (0.d0,0.0d0)
    cpop11 = (0.d0,0.0d0)
    cpop12 = (0.d0,0.0d0)
    cpop21 = (0.d0,0.0d0)
    cpop22 = (0.d0,0.0d0)



    ! Now to read in the basis sets themselves.
    call allocbs(bs1,nbf)
    call readcontbasis(bs1,i1,i2,nbf)
    

    call allocbs(bs2,nbf)
    call readcontbasis(bs2,i3,i4, nbf)
    
  
    !Calculating the overlap between the basis sets and summing over amplitudes
    
    ovrlp12 = ovrlpmatcross(bs1,bs2)
    ovrlp21=transpose(conjg(ovrlp12)) 

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

    crossterm1 = sqrt(nw1*nw2)*real(cpop11+cpop12) 
    crossterm2 = sqrt(nw1*nw2)*real(cpop21+cpop22)

    ! write(6,*) 'pops'
    ! print*,cpop11 
    ! print*,cpop12
    ! print*,cpop21 
    ! print*,cpop22 


    ! write(6,*) 'crossterms, ', crossterm1, crossterm2
   
    
  

  end subroutine cross_terms

  function ovrlpmatcross(bs1,bs2)

    implicit none
    type(basisfn), dimension (:), intent(in)::bs1, bs2
    complex(kind=8), dimension(size(bs1),size(bs1))::output,ovrlpmatcross
    integer::k,j,n

    if (errorflag .ne. 0) return

    ovrlpmatcross = (0.0d0, 0.0d0)

    n = size(bs1)

    do k=1,n
      do j=1,n 
        output(j,k)=ovrlpij(bs2(j)%z,bs1(k)%z)
      end do
    end do
    
    ovrlpmatcross=output
    
    return

  end function ovrlpmatcross

End module clonecondense