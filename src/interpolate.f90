Program Interpolate

! This is a program which interpolates the output data of a set of repeated simulations
! from the MCE program with adaptive stepsizes and also averages the data from 
! the files for each execution subfolder. In its unmodified state this program will
! only work with output created by the MCE program written by CS. 
! The current version number of the MCE program is v0.62 as of 05/02/14

  use neville

  implicit none
  integer :: fileunin, fileunout, cols, filenum, i, j, k, timesteps,ierr, lines, tot, errorflag
  integer, dimension (:), allocatable :: valid
  real(kind=8) :: dt, timeend, timestrt
  real(kind=8), dimension (:), allocatable :: fileline
  real(kind=8), dimension (:,:), allocatable :: input, output, output2
  character (LEN=100) :: filename1, filename2, LINE, LINE1, LINE2, LINE3, LINE4, LINE5, LINE6
  character (LEN=3) :: repstr, folstr
  character (LEN=14) :: myfmt
  logical :: file_exists

  call getarg(0,filename1)
  call getarg(1,LINE)
  if (ierr==-1) then
     print "(a,a)", "Error! Could not read second argument for ", trim(filename1)
     print "(a)", "This should be the number of files in this folder."
     stop
  end if
  read (LINE,*) filenum
  call getarg(2,LINE)
  if (ierr==-1) then
     print "(a,a)", "Error! Could not read third argument for ", trim(filename1)
     print "(a)", "This should be the number of columns in the output files."
     stop
  end if
  read (LINE,*) cols

  OPEN(UNIT=135,file='rundata.csv',status='old',iostat=ierr)

  if (ierr .ne. 0) then
   write(0,"(a)") 'Error in opening rundata.csv file'
   errorflag = 1
   return
  end if

  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)
  read(135,*,iostat=ierr)LINE1, LINE2, LINE3, LINE4, LINE5, LINE6
  if (ierr .ne. 0) then
    write(0,"(a)") 'Error in reading propagation data'
    errorflag = 1
    return
  end if

  close(135)

  read(LINE3,*,iostat=ierr)dt
  if(ierr.ne.0) then
    print *,  "Error reading initial dt"
    stop
  end if
  read(LINE4,*,iostat=ierr)timeend
  if(ierr.ne.0) then
    print *,  "Error reading end time for propagation"
    stop
  end if  
  read(LINE5,*,iostat=ierr)timestrt
  if(ierr.ne.0) then
    print *,  "Error reading starting time for propagation"
    stop
  end if
    
  if ((LINE6(1:1).eq.'S').or.(LINE6(1:1).eq.'s')) then
    print *, "Error! Static step size chosen but interpolation program was started"
    print *, "Exiting now"
    stop
  else if ((LINE6(1:1).ne.'A').or.(LINE6(1:1).ne.'a')) then
    print *, "Error reading step type. Expected 'Static' or 'Adaptive', but got ", trim(LINE2)
    stop
  end if
    
  allocate (valid(filenum))
  valid = 0
  lines = 0
  if (mod(timeend-timestrt,dt).ne.0) timeend=timeend - (mod(timeend-timestrt,dt)) + dt
  timesteps = (timeend-timestrt)/dt + 1
  allocate (output(cols,timesteps))
  output = 0.0d0
  tot = 0

  allocate(fileline(cols))

  do i=1,filenum
     fileunin = 1150 + i
     fileunout = 7150 + i

     lines=0
     
     write(repstr,"(i3.3)") i  

     filename1 = "normpop-"//trim(repstr)//".out"
     filename2 = "normpop-"//trim(repstr)//"_interp.out"

     inquire(file=trim(filename1), exist=file_exists)
     if (file_exists) then
        valid(i) = 1
     else
        valid(i) = 0
     end if

     if (valid(i)==0) cycle
     tot=tot+1

     OPEN(UNIT=fileunin, FILE=trim(filename1),STATUS='OLD', iostat=ierr)
     if (ierr .ne. 0) then
        print *, "Error in opening ", trim(filename1), " file"
        stop
     end if

     OPEN(UNIT=fileunout, FILE=trim(filename2),STATUS='UNKNOWN', iostat=ierr)
     if (ierr .ne. 0) then
        print *, "Error in opening ", trim(filename2), " file"
        stop
     end if

     if (i==1) then
        read(fileunin,"(a141)",iostat=ierr) LINE2
        if (ierr .ne. 0) then
           print *, "Error in reading from ", trim(filename1), " file to get headers"
           stop
        end if
        rewind (fileunin)
     end if

     do while (LINE/="0.00000000E+00")
        read(fileunin,*,iostat=ierr) LINE
        if (ierr .ne. 0) then
           print *, "Error in reading from ", trim(filename1), " file to skip to data"
           stop
        end if
     end do
     backspace (fileunin)
   
     do while (ierr==0)
        LINE3=LINE
        lines = lines + 1
        read(fileunin,*,iostat=ierr) LINE
        if (LINE==LINE3) lines = lines-1
     end do

     lines = lines - 1
     
     allocate(input(cols,lines), stat=ierr)
     if (ierr/=0) then
        print *, "Allocation error in input array"
        stop
     end if
     
     rewind (fileunin)
     do while (LINE/="0.00000000E+000")
        read(fileunin,*,iostat=ierr) LINE
        if (ierr .ne. 0) then
           print *, "Error in reading from ", trim(filename1), " file to skip to data again"
           stop
        end if
     end do
     backspace (fileunin)
   
     ierr=0

     do j=1,lines
        read(fileunin,*,iostat=ierr) fileline(1:cols)
        if (ierr/=0) then
           print *, "Error writing to input array"
           stop
        end if
        do k=1,cols
           input(k,j) = fileline(k)
        end do
     end do

     call nevillealg(input,output,dt,timeend,timestrt,cols,lines,errorflag)   ! Interpolate subroutine for single file
     if (errorflag==1) stop

     write(myfmt,'(a,i0,a)') '(', cols, '(1x,e15.8))'

     do k=1,timesteps
        write(fileunout,myfmt) output(:,k)
     end do     

     rewind (fileunout)

     close (fileunin)

     deallocate(input)

  end do

  deallocate(output)

  OPEN(UNIT=501, FILE="normpop.out",STATUS='new', iostat=ierr) 

  write (501,*), LINE2
  write (501,*), ""
  write (501,*), ""

  allocate (output2(cols,timesteps))

  output2 = 0.0d0
  fileline = 0.0d0

  do i=1,filenum

    if (valid(i)==0) cycle

    fileunout=7150+i

    do k=1,timesteps
       read (fileunout,*,iostat=ierr) fileline
       do j=1,cols
          output2(j,k)=output2(j,k)+fileline(j)
       end do
    end do

    close(fileunout)

  end do

  do k=1,timesteps
     do j=1,cols
        output2(j,k) = output2(j,k)/tot
     end do
  end do 

  write(myfmt,'(a,i0,a)') '(', cols, '(1x,e15.8))'

  do k=1,timesteps
     write(501,myfmt) output2(:,k)
  end do

  close (501)

  deallocate(valid,output2,fileline)

  stop

end program Interpolate
