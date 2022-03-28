Program clonecondense 
    use, intrinsic :: iso_fortran_env
    implicit none
    integer::ierr, j=1, k=1, l=1, n, m, tot=0, reps, cols, folreps, totreps, nlines, i, timestep
    real(kind=8), dimension (:,:), allocatable :: pops
    real(kind=8),dimension(:),allocatable::pop1o, pop2o, popsumo, popdiffo, nrmo, rlacfo, popsav
    real(kind=8),dimension(:),allocatable::imacfo, abacfo, ehro, timeo, rlexo, imexo, abexo
    real(kind=8),dimension(:),allocatable::pop1p, pop2p, popsump, popdiffp, nrmp, rlacfp, popsavp
    real(kind=8),dimension(:),allocatable::imacfp, abacfp, ehrp, timep, rlexp, imexp, abexp
    real(kind=8),dimension(:),allocatable::pop1n, pop2n, popsumn, popdiffn, nrmn, rlacfn
    real(kind=8),dimension(:),allocatable::imacfn, abacfn, ehrn, timen, rlexn, imexn, abexn
    real(kind=8)::pop1av, pop2av, popsumav, popdiffav, nrmav, rlacfav, imacfav
    real(kind=8)::abacfav, rlexav, imexav, abexav, ehrav, timeav
    character(LEN=100)::LINE, filenameo, filenamep, filenamet, original, pair, arg3
    character(LEN=19)::myfmt
    character(LEN=6)::repstr, lstr, timestp
    character(LEN=5000)::lngchar, lngchar2
    integer, dimension(:),allocatable::valid, lines
    logical :: file_exists
    character(LEN=18) :: arg11
    character(LEN=22) :: arg10
    character(len=45) :: command
    call getarg(1,original)
    call getarg(2,pair)
    call getarg(3, arg3)
    read(arg3,*) timestep
   
    

    filenameo =  "normpop-000"//trim(original)//".out"
    filenamep =  "normpop-000"//trim(pair)//".out"
    open(unit=4567, file=trim(filenameo), status='old',iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(3a)") 'Error in opening ', trim(filenameo), ' file'
      stop
    end if
    open(unit=4568, file=trim(filenamep), status='old',iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(3a)") 'Error in opening ', trim(filenamep), ' file'
      stop
    end if
    nlines = 0
    DO
       READ(4567,*,iostat=ierr)
       IF (ierr/=0) EXIT
       nlines = nlines + 1
    END DO
    write(6,*) nlines
    tot = nlines-2
   
    allocate(timeo(tot-1))
    allocate(nrmo(tot-1))
    allocate(rlacfo(tot-1))
    allocate(imacfo(tot-1))
    allocate(abacfo(tot-1))
    allocate(rlexo(tot-1))
    allocate(imexo(tot-1))
    allocate(abexo(tot-1))
    allocate(ehro(tot-1))
    allocate(pop1o(tot-1))
    allocate(pop2o(tot-1))
    allocate(popsumo(tot-1))
    allocate(popdiffo(tot-1))

    allocate(timep(tot-1))
    allocate(nrmp(tot-1))
    allocate(rlacfp(tot-1))
    allocate(imacfp(tot-1))
    allocate(abacfp(tot-1))
    allocate(rlexp(tot-1))
    allocate(imexp(tot-1))
    allocate(abexp(tot-1))
    allocate(ehrp(tot-1))
    allocate(pop1p(tot-1))
    allocate(pop2p(tot-1))
    allocate(popsump(tot-1))
    allocate(popdiffp(tot-1))

    
    allocate(nrmn(tot-1))
    allocate(rlacfn(tot-1))
    allocate(imacfn(tot-1))
    allocate(abacfn(tot-1))
    allocate(rlexn(tot-1))
    allocate(imexn(tot-1))
    allocate(abexn(tot-1))
    allocate(ehrn(tot-1))
    allocate(pop1n(tot-1))
    allocate(pop2n(tot-1))
    allocate(popsumn(tot-1))
    allocate(popdiffn(tot-1))
    
    filenamet = "tempnormpop-000"//trim(original)//".out"
    open(unit=8383, file=filenamet,status='replace')
    write (8383,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
    write (8383,*)""
    write (8383,*)""
    rewind(4567)
    do i=1, tot
      read(4567,*,iostat=ierr) timeo(i-1),nrmo(i-1),rlacfo(i-1),imacfo(i-1),abacfo(i-1),rlexo(i-1),&
                        imexo(i-1),abexo(i-1),ehro(i-1),pop1o(i-1),pop2o(i-1),popsumo(i-1),popdiffo(i-1)
    end do 
    
    do i=1, tot
      read(4568,*,iostat=ierr) timep(i-1),nrmp(i-1),rlacfp(i-1),imacfp(i-1),abacfp(i-1),rlexp(i-1),&
                        imexp(i-1),abexp(i-1),ehrp(i-1),pop1p(i-1),pop2p(i-1),popsump(i-1),popdiffp(i-1)
    end do 
     
    do i=1, timestep
      write(8383,("(13(1x,es16.8e3))"), iostat=ierr) timeo(i),nrmo(i),rlacfo(i),imacfo(i),abacfo(i),rlexo(i),&
                             imexo(i),abexo(i),ehro(i),pop1o(i),pop2o(i),popsumo(i),popdiffo(i)
    end do
    j= timestep +1
    write(6,*)
    do i=timestep, tot-2
      
      nrmn(i) = nrmo(j) + nrmp(j)
      rlacfn(i) = (rlacfo(j)+rlacfp(j))/2
      imacfn(i) =  (imacfo(j)+ imacfp(j))/2
      abacfn(i) = (abacfo(j)+ abacfp(j))/2
      rlexn(i) = (rlexo(j)+ rlexp(j))/2
      imexn(i) = (imexo(j)+ imexp(j))/2
      abexn(i) = (abexo(j)+ abexp(j))/2
      ehrn(i) = (ehro(j) + ehrp(j))/2
      pop1n(i) = pop1o(j) + pop2p(j)
      pop2n(i) = pop2o(j) + pop2p(j)
      write(6,*) "popsum of original = ", popsumo(j)
      write(6,*) "popsum of pair = ", popsump(j)
      popsumn(i) = popsumo(j) + popsump(j)
      write(6,*) "new popsum = ", popsumn(i)
      popdiffn(i) = popdiffo(j) + popdiffp(j)
      j=j+1
    end do 
   
    do i=timestep, tot-2
      write(8383, ("(13(1x,es16.8e3))"),iostat=ierr) timeo(i+1),nrmn(i),rlacfn(i),imacfn(i),abacfn(i),rlexn(i),&
                              imexn(i),abexn(i),ehrn(i),pop1n(i),pop2n(i),popsumn(i),popdiffn(i)
    end do 

    
   

    close(4567)
    close(4568, status= 'delete')
   
    arg10 = "tempnormpop-000"//trim(original)//".out"
    arg11 = "normpop-000"//trim(original)//".out"
    command = "cp "//arg10//arg11
    write(6,*) arg10, arg11, command
    call system(command)
    close(8383)
   


! Can copy from subavrg.f90 for a good starting point?





    
end program 