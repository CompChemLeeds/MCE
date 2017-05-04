Program timehist

  implicit none
  integer istat, n, k, n_box_p, boxn, u, reps, reps2, tot, i, l
  real(kind=8), dimension(:), allocatable :: t
  character(LEN=100) :: LINE, filename
  real(kind=8) :: cutup, cutdown
  integer, dimension(:), allocatable::p_dist, valid
  real:: boxn_rl
  character(LEN=6) :: repstr, lstr
  character(LEN=500) :: lngchar, lngchar2
  logical :: file_exists

  call getarg(0,filename)
  call getarg(1,LINE)
  if (istat==-1) then
    write(0,"(a,a)") "Error! Could not read first argument for ", trim(filename)
    write(0,"(a)") "This should be the number of seed folders used."
    stop
  end if
  read (LINE,*) reps
  call getarg(2,LINE)
  if (istat==-1) then
    write(0,"(a,a)") "Error! Could not read second argument for ", trim(filename)
    write(0,"(a)") "This should be the number of repeats used."
    stop
  end if
  read (LINE,*) reps2

  reps2 = reps2/reps

  allocate(valid(reps))

  do i=1,reps
    write (repstr,"(i0)") i
    filename = "timehist_"//trim(repstr)//".out"    
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
      valid(i) = 1
    else
      valid(i) = 0
    end if
  end do
  tot = sum(valid(:))

  write(6,"(2(i0,a))") tot, " files of ", reps, " valid."

  open (unit=1710,file="timesteps.out",status="unknown",iostat=istat)
  if (istat/=0) then
    write(0,"(a)") "Error in opening timesteps file. Exitting histogram program"
    stop
  end if 
  istat=0
  n=1
  do while (istat==0)
    read (1710,"(a)",iostat=istat) LINE
    n = n+1
  end do
  write(6,"(a,i0)") "size of timestep array is ", n
  rewind(1710)
  allocate(t(n), stat = istat)
  if (istat/=0) then
    write(0,"(a)") "Error in timestep array allocation for histogram program"
    stop
  end if
  do k=1,n
    read (1710,"(e12.5)",iostat=istat) t(k)
    if (istat/=0) cycle
  end do     
  close(1710)

  boxn = 200

  boxn_rl = real(boxn)
  allocate(p_dist(boxn), stat = istat)
  if (istat/=0) then
    write(0,"(a)") "Error in p distribution allocation for histogram"
    stop
  end if

  cutup = maxval(t)
  cutdown = 0.0d0

  open(unit=135,file="timehist_all.out",status="unknown",iostat=istat)
  if (istat/=0) then
    write(0,"(a)") "Error in opening output file for histogram"
    stop
  end if

  do u=1,boxn
    p_dist(u) = 0
  end do

  do u=1,n
    n_box_p = 0      
    if ((t(u).ge.(cutdown-((cutup-cutdown)/(2*boxn_rl)))) &
      .and.(t(u).le.(cutup+((cutup-cutdown)/(2*boxn_rl))))) then
      n_box_p=int((((t(u)-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
      if ((n_box_p.gt.boxn).or.(n_box_p.lt.0)) then
        write(0,"(a,i0)")"Error! Invalid Box Calculated. n_box_p = ", n_box_p
      else if (n_box_p.ne.0) then
        p_dist(n_box_p)=p_dist(n_box_p) + 1
      end if
    end if
  end do

  write (135,*) 'Bin Freq'

  do u=1,boxn
    write (135,"(2(1x,e12.5))") ((cutup-cutdown)*real(u)/boxn_rl)+cutdown, &
                real(p_dist(u))
  end do

  close(135)

  deallocate(t, stat = istat)
  if (istat/=0) then
    write(0,"(a)") "Error in timestep array deallocation for histogram program"
    stop
  end if

  open (unit=176,file="plothistall.gpl",status="unknown",iostat=istat)
  if (istat .ne. 0) then
    write(0,"(a)") 'Error in opening plothistall.gpl output file'
    stop
  end if

  write(176,"(a)") 'set terminal png'
  write(176,"(a)") 'set output "timehist_all.png"'
  write(176,"(a)") 'set title "Histogram of timestep size for all runs combined"'
  write(176,"(a)") 'set xlabel "Timestep"'
  write(176,"(a)") 'set ylabel "Frequency"'
  write(176,"(a)") 'plot "timehist_all.out" u 1:2 t "Combined Histogram" w p'
  close (176)

  if (tot .gt. 1) then

    i=0
    do l=1,reps

      if (valid(l)==0) then
        cycle
      else

        i=i+1
        write(lstr,"(i0)") l
        write(repstr,"(i0)") reps2*i
    
        if (i==1) then     
          lngchar = 'plot "timehist_'//trim(lstr)//'.out" u 1:2 t "'//trim(repstr)//' Reps" w l'
        else
          lngchar2 = lngchar
          lngchar = trim(lngchar2)//', "timehist_'//trim(lstr)//'.out" u 1:2 t "'//trim(repstr)//' Reps" w l'
        end if

      end if

    end do

    open (unit=175,file="plothist.gpl",status="unknown",iostat=istat)
    if (istat .ne. 0) then
      write(0,"(a)") 'Error in opening plothist.gpl output file'
      stop
    end if

    write(175,"(a)") 'set terminal png'
    write(175,"(a)") 'set output "timehist.png"'
    write(175,"(a)") 'set title "Histogram of timestep size for all runs combined"'
    write(175,"(a)") 'set xlabel "Timestep"'
    write(175,"(a)") 'set ylabel "Frequency"'
    write(175,"(a)") trim(lngchar)
    close (175)

  end if

  stop

end program timehist
