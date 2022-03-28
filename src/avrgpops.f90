Program avrgnorm

  ! This program averages a set of normpop.out files from different instances of the
  ! program, outputting cumulativly averaged files (used to check for convergence)
  ! and also residuals. The program also creates a set of gnuplot files which are
  ! used to plot the combined results. This program should be called from the 
  ! <<collate.sh>> script, and should be used with care if called manually.
  
  ! The arguments for the program are:
  !     1: The number of seed folders
  !     2: The total number of repeats
  !     3: The total number of columns in the normpop.out files

  implicit none

  real(kind=8), dimension(:,:), allocatable :: pops
  real(kind=8), dimension(:), allocatable :: pop1, pop2, popsum, popdiff, nrm, rlacf 
  real(kind=8), dimension(:), allocatable :: popsav, imacf, abacf, ehr, time, rlex
  real(kind=8), dimension(:), allocatable :: imex, abex
  real(kind=8)::pop1av, pop2av, popsumav, popdiffav, nrmav, rlacfav, imacfav
  real(kind=8)::abacfav, rlexav, imexav, abexav, ehrav, timeav
  integer, dimension(:), allocatable :: valid, lines
  integer::ierr, i=1, j=1, k=1, l=1, n, m, tot=0, folders, folreps, totreps, cols, v
  character(LEN=50000)::lngchar, lngchar2
  character(LEN=100)::LINE, filename
  character(LEN=19)::myfmt
  character(LEN=6)::repstr, lstr
  logical :: file_exists
  
  ! Read in the preprocessor arguments
  call getarg(0,filename)
  call getarg(1,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read first argument for ", trim(filename)
    write (0,"(a)") "This should be the number of seed folders used."
    stop
  end if
  read (LINE,*)folders
  call getarg(2,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read second argument for ", trim(filename)
    write (0,"(a)") "This should be the number of repeats used."
    stop
  end if
  read (LINE,*) totreps
  call getarg(3,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read third argument for ", trim(filename)
    write (0,"(a)") "This should be the number of columns in each normpop file."
    stop
  end if
  read (LINE,*) cols  

  folreps = totreps/folders
  allocate(valid(folders))
  allocate(lines(folders))
  lines = 0

  ! Check to see which folders have returned valid normpop.out files
  ! assigning array values for later reference
  do i=1,folders
    write (repstr,"(i0)") i
    filename = "normpop_"//trim(repstr)//".out"    
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
      valid(i) = 1
    else
      valid(i) = 0
    end if
  end do

  ! Open the normpop.out files and count the number of lines each has
  do i=1,folders       
    if (valid(i)==0) then
      cycle
    end if
    tot = tot+1       
    write (repstr,"(i0)") i
    filename = "normpop_"//trim(repstr)//".out"
    OPEN(UNIT=128, FILE=trim(filename),STATUS='OLD', iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(3a,i0)") 'Error in opening ',trim(filename),' file. Ierr was ', ierr
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
  
  ! Check that each normpop.out file has the same number of lines
  ! As any interpolation will have already been carried out at this stage, this
  ! always be the case for successful executions. Mismatches could indicate that
  ! the simulations are incomplete and need restarting
  do i=1,folders
    if ((lines(i).ne.lines(1)).and.(lines(i).ne.0)) then
      write (0,"(a,i0)") "Error - Line number mismatch in repeat ", i
      write (0,"(a,i0,a,i0)") "Expected ", lines(1), " lines, but got ", lines(i) 
      stop
    end if
  end do

  write (6,"(i0,a,i0,a)") tot, " files of ",folders, " were valid."
  
  ! Allocate the arrays for each column. If cols=13 then npes=2, so the population
  ! sum and population difference will be output also
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

  n = 1130  !For file units
  m = 7150

  ! Open each file again, and advance read pointer to beginning of data
  do i=1,folders
       
    if (valid(i)==0) then
      cycle
    end if     
    write (repstr,"(i0)") i
    filename = "normpop_"//trim(repstr)//".out"
    n = n+1
    OPEN(UNIT=n, FILE=trim(filename),STATUS='OLD', iostat=ierr)

    if (ierr .ne. 0) then
      write (0,"(3a,i0)") 'Error in opening ',trim(filename),' file. Ierr was ', ierr
      stop
    end if
    v=1
    do
      read (n,*,iostat=ierr)LINE
      write(6,*) v
      v = v+1
      if (ierr.ne.0) then

        write (0,"(3a,i0)") "Read error in normpop_",trim(repstr),".out. Ierr was ",&
                            ierr
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
  ! Open the output files and write the headers
  do l=1,tot

    m=7150+l

    write(repstr,"(i0)") folreps*l

    filename= "normpop_cumul_"//trim(repstr)//".out"

    open (unit=m, file=trim(filename), status='unknown', iostat=ierr)

    if (ierr .ne. 0) then
      write (0,"(4a,i0)") 'Error in opening ', trim(filename), ' output file.', &
                        ' Ierr was ', ierr
      stop
    end if

    if (cols==13) then
      write (m,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) ", &
                          "Im(Extra) |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
    else
      write (m,"(2a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) ", &
                          "Im(Extra) |Extra| Sum(HEhr) Pops(1..n)"
    end if
    write (m,*) ""
    write (m,*) ""

  end do


  do k = 1,lines(1)

    ! Read a single line from each of the normpop.out files
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

    ! Calculate the cumulative averages and write them to file
    do i=1,tot

      m=7150+i
       
      timeav = sum(time(1:i))/i
      nrmav = sum(nrm(1:i))/i
      rlacfav = sum(rlacf(1:i))/i
      imacfav = sum(imacf(1:i))/i
      abacfav = sum(abacf(1:i))/i
      rlexav = sum(rlex(1:i))/i
      imexav = sum(imex(1:i))/i
      abexav = sum(abex(1:i))/i
      ehrav = sum(ehr(1:i))/i
      if (cols==13) then
        pop1av = sum(pop1(1:i))/i
        pop2av = sum(pop2(1:i))/i
        popsumav = sum(popsum(1:i))/i
        popdiffav = sum(popdiff(1:i))/i
      else
        do j=1,size(popsav)
          popsav(j) = sum(pops(j,1:i))/i
        end do
      end if

      if (cols == 13) then
        write (m,"(13(1x,es16.8e3))") timeav,nrmav,rlacfav,imacfav,abacfav,rlexav,&
                        imexav,abexav,ehrav,pop1av,pop2av,popsumav,popdiffav
      else
        write(myfmt,'( "(",i0,"(1x,es16.8e3))" )') cols
        write(m,myfmt) timeav,nrmav,rlacfav,imacfav,abacfav,rlexav,&
                        imexav,abexav,ehrav,popsav
      end if

    end do

  end do

  ! Create the gnuplot script for comparison of the population differences
  do l=1,tot

    write(repstr,"(i0)") folreps*l
    
    if (l==1) then     
      lngchar = 'plot "normpop_cumul_'//trim(repstr)//'.out" u 1:13 t "' & 
                            //trim(repstr)//' Reps" w l'
    else
      lngchar2 = lngchar
      lngchar = trim(lngchar2)//', "normpop_cumul_'//trim(repstr)// & 
                       '.out" u 1:13 t "'//trim(repstr)//' Reps" w l'
    end if

  end do

  lngchar2 = lngchar
  lngchar = trim(lngchar2)//', "normpop_cumul_'//trim(repstr)// & 
                                      '.out" u 1:2 t "Total Av Norm" w l'

  open (unit=175,file="plotpopdiff.gpl",status="unknown",iostat=ierr)
  if (ierr .ne. 0) then
    write (0,"(2a,i0)") 'Error in opening plotpopdiff.gpl output file.', &
                              ' Ierr was ', ierr
    stop
  end if

  write(175,"(a)") 'set terminal png'
  write(175,"(a)") 'set output "popsculm.png"'
  write(175,"(a)") 'set title "Graph of convergence over multiple repetitions"'
  write(175,"(a)") 'set xlabel "Time"'
  write(175,"(a)") 'set ylabel "Population difference"'
  write(175,"(a)") trim(lngchar)
  close (175)
  
  ! Create gnuplot script for total averaged population difference
  open (unit=175,file="plottotpopdiff.gpl",status="unknown",iostat=ierr)
  if (ierr .ne. 0) then
    write (0,"(2a,i0)") 'Error in opening plottotpopdiff.gpl output file.', &
                              ' Ierr was ', ierr
    stop
  end if

  write(175,"(a)") 'set terminal png'
  write(175,"(a)") 'set output "popsdifftot.png"'
  write(175,"(a)") 'set title "Graph total population difference"'
  write(175,"(a)") 'set xlabel "Time"'
  write(175,"(a)") 'set ylabel "Population difference"'
  write(175,"(a)") 'plot "normpop_cumul_'//trim(repstr)//'.out" u 1:13 t "' & 
                         //trim(repstr)//' Reps" w l, "" u 1:2 t "Total Av Norm" w l'
  close (175)  
  
  ! Create gnuplot script for total averaged populations
  open (unit=175,file="plottotpops.gpl",status="unknown",iostat=ierr)
  if (ierr .ne. 0) then
    write (0,"(2a,i0)") 'Error in opening plottotpops.gpl output file.', &
                              ' Ierr was ', ierr
    stop
  end if

  write(175,"(a)") 'set terminal png'
  write(175,"(a)") 'set output "popstot.png"'
  write(175,"(a)") 'set title "Graph total populations"'
  write(175,"(a)") 'set xlabel "Time"'
  write(175,"(a)") 'set ylabel "Populations"'
  write(175,"(a)") 'plot "normpop_cumul_'//trim(repstr)//'.out" u 1:10 t "' & 
                         //trim(repstr)//' Reps" w l, "" u 1:11 t "'//trim(repstr)//'" w l'
  close (175) 
  
  ! Create gnuplot script for extra calculated quantity (ie dispersion or dipole acceleration)
  open (unit=175,file="plotext.gpl",status="unknown",iostat=ierr)
  if (ierr .ne. 0) then
    write (0,"(2a,i0)") 'Error in opening plotext.gpl output file.', &
                              ' Ierr was ', ierr
    stop
  end if

  write(175,"(a)") 'set terminal png'
  write(175,"(a)") 'set output "ext.png"'
  write(175,"(a)") 'set title "Graph total extra quantity"'
  write(175,"(a)") 'set xlabel "Time"'
  write(175,"(a)") 'set ylabel "Extra Quantity"'
  write(175,"(a)") 'plot "normpop_cumul_'//trim(repstr)//'.out" u 1:6 t "' & 
                         //trim(repstr)//' Reps" w l'
  close (175)

  ! Create the population difference residuals file
  if ((tot .gt. 1).and.(cols==13)) then

    open(unit=180,file="popdiffresiduals.out",status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(2a,i0)") 'Error in opening popdiffresiduals.out output file.', &
                              ' Ierr was ', ierr
      stop
    end if

    do i=1,tot
      n=7150+i
      rewind(n)
      do
        read(n,*,iostat=ierr) LINE
        if (LINE=="0.00000000E+000") then
          backspace(n)
          exit
        else
          cycle
        end if
      end do          
    end do

    do k = 1,lines(1)

      ! Read in data from the cumulative average files
      do i=1,tot
        n=7150+i
        read(n,*,iostat=ierr)time(i),nrm(i),rlacf(i),imacf(i),abacf(i),rlex(i),&
                      imex(i),abex(i),ehr(i),pop1(i),pop2(i),popsum(i),popdiff(i)
        if (time(i).ne.time(1)) then
          write (0,"(a,a,es16.8e3,a,i0)"),"File synchronisation error when reading",&
                               " at time ", time(1) ," in unit ", n
          write (0,"(a,i0,a,i0,a,i0)") "i = ", i, " n = ", n
          stop
        end if
      end do

      ! Calculate residual and write to file
      do i=1,tot
        popdiff(i)=popdiff(tot)-popdiff(i)
      end do
      write (myfmt,'(a,i4.4,a)') '(', tot+1, '(1x,es16.8e3))'
      write (180,trim(myfmt)) time(1), popdiff

    end do

    close (180)

    ! Create the gnuplot script for plotting the population difference residuals
    do l=1,tot

      write(repstr,"(i0)") folreps*l
      write(lstr,"(i0)") l+1
    
      if (l==1) then     
        lngchar = 'plot "popdiffresiduals.out" u 1:2 t "'//trim(repstr)//&
                     ' Reps" w l'
      else
        lngchar2 = lngchar
        lngchar = trim(lngchar2)//', "" u 1:'//trim(lstr)//' t "'//trim(repstr)//&
                     ' Reps" w l'
      end if

    end do

    open (unit=176,file="plotpopres.gpl",status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(a,i0)") 'Error in opening plotpopres.gpl output file. Ierr was',ierr
      stop
    end if

    write(176,"(a)") 'set terminal png'
    write(176,"(a)") 'set output "popsculmres.png"'
    write(176,"(a)") 'set title "Graph of residuals over multiple repetitions"'
    write(176,"(a)") 'set xlabel "Time"'
    write(176,"(a)") 'set ylabel "Residual Population Difference"'
    write(176,"(a)") trim(lngchar)
    close (176)
  end if
  
  write (6,*) "Averages completed and gpl files created" 

  stop

end program avrgnorm
