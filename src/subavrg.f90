Program subavrg

  implicit none
  integer::ierr, i=1, j=1, k=1, l=1, n, m, tot=0, reps, cols, folreps, totreps
  real(kind=8), dimension (:,:), allocatable :: pops
  real(kind=8),dimension(:),allocatable::pop1, pop2, popsum, popdiff, nrm, rlacf, popsav
  real(kind=8),dimension(:),allocatable::imacf, abacf, ehr, time, rlex, imex, abex
  real(kind=8)::pop1av, pop2av, popsumav, popdiffav, nrmav, rlacfav, imacfav
  real(kind=8)::abacfav, rlexav, imexav, abexav, ehrav, timeav
  character(LEN=100)::LINE, filename
  character(LEN=19)::myfmt
  character(LEN=6)::repstr, lstr, timestp
  character(LEN=5000)::lngchar, lngchar2
  integer, dimension(:),allocatable::valid, lines
  logical :: file_exists
      
  call getarg(0,filename)
  call getarg(1,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read first argument for ", trim(filename)
    write (0,"(a)") "This should be the number of repeats used."
    stop
  end if
  read (LINE,*) reps
  call getarg(2,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read second argument for ", trim(filename)
    write (0,"(a)") "This should be the number of columns per file used."
    stop
  end if
  read (LINE,*) cols

  allocate(valid(reps))
  allocate(lines(reps))
  lines = 0

  do i=1,reps
    write (repstr,"(i4.4)") i
    filename = "normpop-"//trim(repstr)//".out"    
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
      valid(i) = 1
      print*,1
    else
      valid(i) = 0
    end if
  end do

  do i=1,reps       
    if (valid(i)==0) then
      cycle
    end if
    tot = tot+1       
    write (repstr,"(i4.4)") i
    filename = "normpop-"//trim(repstr)//".out"
    OPEN(UNIT=128, FILE=trim(filename),STATUS='OLD', iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(3a)") 'Error in opening ', trim(filename), ' file'
      stop
    end if
    read(128,*,iostat=ierr) LINE
    do while (ierr==0)
      lines(i) = lines(i) + 1
      read(128,*,iostat=ierr) LINE
    end do
    close (128)
    lines(i) = lines(i) - 1
  end do
 
  call getcwd(LINE)
  do i=1,reps
    if ((lines(i).ne.lines(1)).and.(lines(i).ne.0)) then
      write (0,"(a,i0)") "Error - Line number mismatch in repeat ", i
      write (0,"(a,i0,a,i0)") "Expected ", lines(1), " lines, but got ", lines(i)
      write (0,"(2a)") "This occured in ", trim(LINE)
      stop
    end if
  end do
       
  allocate(time(tot))
  allocate(nrm(tot))
  allocate(rlacf(tot))
  allocate(imacf(tot))
  allocate(abacf(tot))
  allocate(rlex(tot))
  allocate(imex(tot))
  allocate(abex(tot))
  allocate(ehr(tot))
  if (cols==13) then
    allocate(pop1(tot))
    allocate(pop2(tot))
    allocate(popsum(tot))
    allocate(popdiff(tot))
  else
    allocate(pops(cols-9,tot))
    allocate(popsav(cols-9))
  end if
  n = 1130
  m = 7150

  do i=1,reps
       
    if (valid(i)==0) then
      cycle
    end if     
    write (repstr,"(i4.4)") i
    filename = "normpop-"//trim(repstr)//".out"
    n = n+1
    OPEN(UNIT=n, FILE=trim(filename),STATUS='OLD', iostat=ierr)

    if (ierr .ne. 0) then
      write (0,"(a,a,a)") 'Error in opening ', trim(filename), ' file'
      stop
    end if
    do
      read (n,*,iostat=ierr)LINE
      if (ierr.ne.0) then
        write (0,"(a,a,a)") "Read Error in normpop-", trim(repstr), ".out"
        stop
      end if
      if (LINE=="0.00000000E+000") then
        backspace (n)
        exit
      else
        cycle
      end if
    end do

  end do

  filename="normpop.out"

  open (unit=m, file=trim(filename), status='unknown', iostat=ierr)

  if (ierr .ne. 0) then
    write (0,"(a)") 'Error in opening normpop.out output file'
    stop
  end if

  if (cols==13) then
    write (m,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
  else
    write (m,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pops(1..n)"
  end if
  write (m,*)""
  write (m,*)""

  do k = 1,lines(1)

    do i=1,tot
      n=1130+i
      if (cols==13) then
        read(n,*,iostat=ierr)time(i),nrm(i),rlacf(i),imacf(i),abacf(i),rlex(i),&
                        imex(i),abex(i),ehr(i),pop1(i),pop2(i),popsum(i),popdiff(i)
      else
        read(n,*,iostat=ierr)time(i),nrm(i),rlacf(i),imacf(i),abacf(i),rlex(i),&
                        imex(i),abex(i),ehr(i),pops(:,i)
      end if     
      if (time(i).ne.time(1)) then
        write (0,"(a,e15.8,a,i0)") "File synchronisation error when reading at ", & 
                                    "time ", time(1) ," in unit ", n
        write (0,"(a,i0,a,i0,a,i0)") "i = ", i, " n = ", n
        write (0,"(a,e15.8,a,e15.8)") "Expected t=", time(1), " but got t=", time(i)
        stop
      end if
    end do

    timeav = sum(time(:))/tot
    nrmav = sum(nrm(:))/tot
    rlacfav = sum(rlacf(:))/tot
    imacfav = sum(imacf(:))/tot
    abacfav = sum(abacf(:))/tot
    rlexav = sum(rlex(:))/tot
    imexav = sum(imex(:))/tot
    abexav = sum(abex(:))/tot
    ehrav = sum(ehr(:))/tot
    if (cols==13) then
      pop1av = sum(pop1(:))/tot
      pop2av = sum(pop2(:))/tot
      popsumav = sum(popsum(:))/tot
      popdiffav = sum(popdiff(:))/tot
    else
      do j=1,size(popsav)
        popsav(j) = sum(pops(j,:))
      end do
    end if
    if (cols == 13) then
      write (m,"(13(1x,es16.8e3))") timeav,nrmav,rlacfav,imacfav,abacfav,rlexav,&
                      imexav,abexav,ehrav,pop1av,pop2av,popsumav,popdiffav
    else
      write(myfmt,"(a,i0,a)") '"(', cols, '(1x,es16.8e3))"'
      write(m,myfmt) timeav,nrmav,rlacfav,imacfav,abacfav,rlexav,&
                      imexav,abexav,ehrav,popsav
    end if
  end do
  
  close (m)

  stop

end program subavrg
