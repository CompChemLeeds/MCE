MODULE outputs

   use globvars

!*************************************************************************************************!
!*
!*         Output Module
!*           
!*   Contains subroutines for:
!*
!*      1) Generating and outputting a histogram using a complex variable $$
!*      2) Generating and outputting a histogram using a real variable
!*      3) Generating and outputting a two-dimensional histogram using a complex variable $$
!*      4) Outputting the basis set values and parameters for later use
!*      5) Creating a file with headers for outputting the simulation data for a single thread 
!*      6) Outputting the simulation data at a single timestep for a single thread 
!*      7) Outputting the simulation data for all timesteps over all threads
!*      8) Creating a file with headers and plotting scripts for outputting the wavefunction values for a single thread
!*      9) Outputting the wavefunction values at a single timestep for a single thread
!*     10) Creating a file with headers and plotting scripts for outputting the trajectories for a single thread
!*     11) Outputting the trajectory values at a single timestep for a single thread
!*     12) Interpolation controller subroutine for averaging all adaptive timestep data for a single multi-thread run
!*     13) Neville Algorithm subroutine which controls the polynomial interpolation for a single thread
!*     14) Polynomial Interpolation subroutine which calculates the value of a single datapoint **
!*     15) Data seeking algorithm which finds the array index of a particular time value in an ordered list **
!*     
!*      $$ - this subroutine is not currently used and was created for debugging purposes
!*      ** - this subroutine was copied verbatim from Press et al Numerical Recipes.
!*      
!*************************************************************************************************!

contains

!*************************************************************************************************!
!        Output Subroutines
!*************************************************************************************************!

   subroutine histogram(A, n, filename, cutup, cutdown) ! Level 1 Subroutine

    implicit none

    integer:: n_box_p, boxn, u, ierr
    integer, dimension(:), allocatable::p_dist
    real:: boxn_rl
    real(kind=8),intent(inout)::cutup, cutdown
    real(kind=8),dimension(:),intent(in)::A
    integer,intent(in)::n
    character(LEN=12),intent(in)::filename

    if (errorflag .ne. 0) return

    ierr = 0
    boxn = 200

    boxn_rl = real(boxn)
    allocate(p_dist(boxn), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in p distribution allocation for histogram"
      errorflag=1
      return
    end if

    if ((cutup == 0.0).and.(cutdown == 0.0)) then
      cutup = maxval(A)
      cutdown = minval(A)
    end if

    open(unit=135,file=filename,status='unknown',iostat=ierr)

    do u=1,boxn
      p_dist(u) = 0
    end do

    do u=1,n
      n_box_p = 0      
      if ((A(u).ge.(cutdown-((cutup-cutdown)/(2*boxn_rl)))) &
         .and.(A(u).le.(cutup+((cutup-cutdown)/(2*boxn_rl))))) then
        n_box_p=int((((A(u)-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
        if ((n_box_p.gt.boxn).or.(n_box_p.lt.0)) then
          write(0,"(a)")"Error! Invalid Box Calculated. n_box_p = ", n_box_p
        else if (n_box_p.ne.0) then
          p_dist(n_box_p)=p_dist(n_box_p) + 1
        end if
      end if
    end do

    write (135,*) 'Bin Freq'

    do u=1,boxn
      write (135,"(2(1x,e13.5e3))") ((cutup-cutdown)*real(u)/boxn_rl)+cutdown, &
                  real(p_dist(u))
    end do

    close(135)

    return

   end subroutine histogram

!--------------------------------------------------------------------------------------------------

   subroutine outbs(bs, reps, mup, muq, t, x)   !   Level 1 Subroutine
    implicit none
    type(basisfn), dimension (:), intent(in) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: t
    integer::m, j, r, k, ierr, bsunit, fileun
    character(LEN=22)::filename, filenm, filenm2, myfmt
    integer, intent(in) :: reps, x
    character(LEN=4):: rep
    character(LEN=5)::step
    character(LEN=1)::rkstp

    ierr = 0

    write(rep,"(i4.4)") reps
    if (x.lt.0) then
      write(step,"(i5.4)") x
    else
      write(step,"(i5.5)") x
    end if
    
    if (errorflag.eq.0) then
      filename = "Outbs-"//trim(rep)//".out"
    else
      filename = "Errbs-"//trim(rep)//".out" 
    end if

    bsunit=135+reps

    open(unit=bsunit, file=trim(filename), status='unknown', iostat=ierr)
    
    if (size(bs).eq.0) then
      write(bsunit,"(a)") "Basis set cannot be output as it is not allocated properly"
      errorflag = 1
      return
    end if

    write(bsunit,"(a,1x,i4)"       ) 'ndof'       , ndim
    write(bsunit,"(a,1x,i4)"       ) 'nconf'      , npes
    write(bsunit,"(a,1x,i4)"       ) 'nbasisfns'  , size(bs)
    write(bsunit,"(a,1x,i4)"       ) 'initial_PES', in_pes
    write(bsunit,"(a,1x,es25.17e3)") 'time'       , t
    write(bsunit,*),""
    do m=1,ndim
      write(bsunit,"(a,1x,i3,2(1x,es25.17e3))") 'zinit', m, muq(m), mup(m)
    end do

    do j=1,size(bs)
      write(bsunit,*),""
      write(bsunit,"(a,1x,i5)")'basis', j
      write(bsunit,"(a,2(1x,es25.17e3))") "D", dble(bs(j)%D_big), dimag(bs(j)%D_big)
      do r=1,npes
        write(bsunit,"(a,i4,2(1x,es25.17e3))") "a",r, dble(bs(j)%a_pes(r)), dimag(bs(j)%a_pes(r)) 
      end do  
      do r=1,npes
        write(bsunit,"(a,i4,2(1x,es25.17e3))") "d",r, dble(bs(j)%d_pes(r)), dimag(bs(j)%d_pes(r)) 
      end do  
      do r=1,npes
        write(bsunit,"(a,i4,1x,es25.17e3)") "s",r, (bs(j)%s_pes(r))
      end do  
      do m=1,ndim
        write(bsunit,"(a,i4,2(1x,es25.17e3))") "z",m, dble(bs(j)%z(m)), dimag(bs(j)%z(m))
      end do
    end do

    close(bsunit)
   
    return

   end subroutine outbs

!--------------------------------------------------------------------------------------------------

  subroutine outclones(clonenum, reps, clone)
  
    implicit none
    
    integer, dimension(:), intent(in) :: clonenum, clone
    integer, intent(in) :: reps
    
    integer :: j, n, ierr
    character(LEN=255) :: filename
    character(LEN=4) :: rep
    
    
    if (errorflag .ne. 0) return
    
    write(rep,"(i4.4)") reps
    filename="clonearr-"//rep//".out"
    
    n = 90254+reps
    
    open(unit=n,file=trim(filename),status="unknown",iostat=ierr)
    if (ierr/=0) then
      write (0,"(3a,i0)") "Error opening ", trim(filename), " file. Ierr was ", ierr
      errorflag = 1
      return
    end if

    do j=1,size(clonenum)
      write(n,"(i5,1x,i2,1x,i5)") j, clonenum(j), clone(j)
    end do
    
    close(n)
    
  end subroutine outclones

!--------------------------------------------------------------------------------------------------

   subroutine outnormpopadapheads(reps)   !   Level 1 Subroutine

    implicit none
    integer :: ierr, fileun
    character(LEN=16)::filenm, filenm2
    integer, intent(in) :: reps
    character(LEN=4):: rep
    logical :: file_exists

    if (errorflag .ne. 0) return

    ierr = 0

    write(rep,"(i4.4)") reps

    filenm = "normpop-"//trim(rep)//".out"  
    inquire(file=filenm,exist=file_exists)
    if(file_exists.eqv..false.) then
      fileun = 610+reps

      open(unit=fileun,file=filenm,status='unknown',iostat=ierr)

      if (ierr.ne.0) then
        write(0,"(a,a,a)") "Error in initial opening of ", trim(filenm), " file"
        errorflag = 1
        return
      end if

      if (npes==2) then 
        write(fileun,"(a,a)") "Time Norm Re(ACF(t)) Im(ACF(t)) |ACF(t)| Re(Extra) Im(Extra) |Extra| ",&
                        "Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
        write(fileun,*), ""
        write(fileun,*), ""
      else
        write(fileun,"(a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pops(1...n)"
        write(fileun,*), ""
        write(fileun,*), ""
      end if

      close(fileun)
    end if

    filenm2="plotacf-"//rep//".gnu"
    inquire(file=filenm2,exist=file_exists)
    if(file_exists.eqv..false.) then
      open (unit=fileun,file=filenm2,status="new",iostat=ierr)
      if (ierr .ne. 0) then
        write(0,"(a,a,a)") 'Error in opening ', trim(filenm2), ' output file'
        errorflag=1
        return
      end if

      write(fileun,"(a)") 'set terminal png'
      write(fileun,"(a,a,a)") 'set output "ACF-',trim(rep),'.png"'
      write(fileun,"(a)") 'set title "Graph of Autocorrelation Function"'
      write(fileun,"(a)") 'set xlabel "Time (au)"'
      write(fileun,"(a)") 'set ylabel "ACF"'
      write(fileun,"(a,a,a)") 'plot "', filenm, '" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l'
      close(fileun)
    end if

    filenm2="plotnrm-"//rep//".gnu"
    inquire(file=filenm2,exist=file_exists)
    if(file_exists.eqv..false.) then
      open (unit=fileun,file=filenm2,status="new",iostat=ierr)
      if (ierr .ne. 0) then
        write(0,"(a,a,a)") 'Error in opening ', trim(filenm2), ' output file'
        errorflag=1
        return
      end if

      write(fileun,"(a)") 'set terminal png'
      write(fileun,"(a,a,a)") 'set output "Norm-',trim(rep),'.png"'
      write(fileun,"(a)") 'set title "Graph of Norm"'
      write(fileun,"(a)") 'set xlabel "Time (au)"'
      write(fileun,"(a)") 'set ylabel "Norm"'
      write(fileun,"(a,a,a)") 'plot "', filenm, '" u 1:2 t "Norm" w l'
      close(fileun)
    end if

    filenm2="plotext-"//rep//".gnu"
    inquire(file=filenm2,exist=file_exists)
    if(file_exists.eqv..false.) then
      open (unit=fileun,file=filenm2,status="unknown",iostat=ierr)
      if (ierr .ne. 0) then
        write(0,"(a,a,a)") 'Error in opening ', trim(filenm2), ' output file'
        errorflag=1
        return
      end if

      write(fileun,"(a)") 'set terminal png'
      if (sys=="IV") then
        write(fileun,"(a,a,a)") 'set output "Dipole-',trim(rep),'.png"'
        write(fileun,"(a)") 'set title "Graph of Dipole Acceleration"'
        write(fileun,"(a)") 'set ylabel "Dipole Acceleration"'
        write(fileun,"(a)") 'set xlabel "Time (au)"'
        write(fileun,"(a,a,a)") 'plot "', filenm, '" u 1:6 t "Real" w l'
      else if (sys=="VP") then
        write(fileun,"(a,a,a)") 'set output "Disp-',trim(rep),'.png"'
        write(fileun,"(a)") 'set title "Graph of Dispersion"'
        write(fileun,"(a)") 'set ylabel "<\Delta x>"'
        write(fileun,"(a)") 'set xlabel "Time (au)"'
        write(fileun,"(a,a,a)") 'plot "', filenm, '" u 1:6 t "Real" w l'
      else
        write(fileun,"(a,a,a)") 'set output "Extra-',trim(rep),'.png"'
        write(fileun,"(a)") 'set title "Graph of Extra Calculated Quantity"'
        write(fileun,"(a)") 'set ylabel "Extra"'
        write(fileun,"(a)") 'set xlabel "Time (au)"'
        write(fileun,"(a,a,a)") 'plot "', filenm, '" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l'
      end if
      
      close(fileun)
    end if


    filenm2="plotdif-"//rep//".gnu"
    inquire(file=filenm2,exist=file_exists)
    if(file_exists.eqv..false.) then
      open (unit=fileun,file=filenm2,status="unknown",iostat=ierr)
      if (ierr .ne. 0) then
        write(0,"(a,a,a)") 'Error in opening ', trim(filenm2), ' output file'
        errorflag=1
        return
      end if

      write(fileun,"(a)") 'set terminal png'
      write(fileun,"(a,a,a)") 'set output "PopDiff-',trim(rep),'.png"'
      write(fileun,"(a)") 'set title "Graph of Population Difference"'
      write(fileun,"(a)") 'set ylabel "Population Difference"'
      write(fileun,"(a)") 'set xlabel "Time (au)"'
      write(fileun,"(a,a,a)") 'plot "', filenm, '" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l'
      close(fileun)
    end if

   end subroutine outnormpopadapheads

!--------------------------------------------------------------------------------------------------

   subroutine outnormpopadap(norm, acft, extra, ehr, popt, x, reps, time)   !   Level 1 Subroutine   

    implicit none
    real(kind=8), intent (in) :: norm, ehr, time
    complex(kind=8), intent (in) :: acft, extra
    real(kind=8), dimension (:), intent (in) :: popt
    integer, intent (in) :: x
    integer :: ierr, fileun
    character(LEN=18) :: filenm,myfmt
    integer, intent(in) :: reps
    character(LEN=4):: rep

    if (errorflag .ne. 0) return

    ierr = 0

    write(rep,"(i4.4)") reps

    filenm = "normpop-"//trim(rep)//".out"     

    fileun = 610+reps  

    open(unit=fileun,file=trim(filenm),status='old',access='append',iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a,a,a,1x,i5)") "Error opening ", trim(filenm), " file in step", x
      errorflag = 1
      return
    end if

    write(myfmt,'(a,i0,a)') '(', 9+npes, '(1x,e16.8e3))'

    if (npes==2) then
      write(fileun,"(13(1x,e16.8e3))") time, norm, dble(acft), dimag(acft), abs(acft), &
                     dble(extra), dimag(extra), abs(extra), ehr, &
                     popt(1), popt(2), popt(1)+popt(2), popt(1)-popt(2)
    else
       write(fileun,myfmt), time, norm, dble(acft), dimag(acft), abs(acft), &
                     dble(extra), dimag(extra), abs(extra), ehr, popt(:)
    end if
    

    close (fileun)
   
   end subroutine outnormpopadap

!--------------------------------------------------------------------------------------------------

   subroutine outnormpopstat(norm, acft, extra, ehr, pops)   !   Level 1 Subroutine

    implicit none
    real(kind=8),dimension(:), intent (in) :: norm, ehr
    complex(kind=8), dimension(:), intent (in) :: acft, extra
    real(kind=8),dimension(:,:), intent(in) :: pops
    real(kind=8) :: time
    integer :: ierr, t, arrend
    character(LEN=11) :: filenm
    character(LEN=17) :: myfmt

    if (errorflag .ne. 0) return 

    ierr = 0   

    filenm = "normpop.out"

    open(unit=150,file=trim(filenm),status='unknown',iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a,1x,i5)") "Error opening ", trim(filenm), " file"
      errorflag = 1
      return
    end if

    if (npes==2) then 
      write(150,"(a,a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra) |Extra| ", &
                    "Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
      write(150,*), ""
      write(150,*), ""
    else
      write(150,"(a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pops(1...n)"
      write(150,*), ""
      write(150,*), ""
      write(myfmt,'(a,i0,a)') '(', 9+npes, '(1x,e16.8e3))'
    end if
    
    arrend = size(norm) - 1
    
    if (npes==2) then
      do t=1,arrend
        time = dtinit * (t-1)
        write(150,"(13(1x,e16.8e3))") time, norm(t), dble(acft(t)), dimag(acft(t)), abs(acft(t)), &
              dble(extra(t)), dimag(extra(t)), abs(extra(t)), ehr(t), &
              pops(t,1), pops(t,2), pops(t,1)+pops(t,2), pops(t,1)-pops(t,2)
      end do
    else
      do t=1,arrend
        time = dtinit * (t-1)
        write(150,myfmt), time, norm(t), dble(acft(t)), dimag(acft(t)), abs(acft(t)), &
              dble(extra(t)), dimag(extra(t)), abs(extra(t)), ehr(t), pops(t,:)
      end do
    end if

    close (150)
   
   end subroutine outnormpopstat

!--------------------------------------------------------------------------------------------------

   subroutine interpolate(cols,errorflag)

    implicit none
    integer, intent (inout) :: errorflag
    integer, intent (in) :: cols
    integer :: fileunin, fileunout, i, j, k, timesteps,ierr, lines, tot
    integer, dimension (:), allocatable :: valid
    real(kind=8), dimension (:), allocatable :: fileline
    real(kind=8), dimension (:,:), allocatable :: input, output, output2
    character (LEN=100) :: filename1, filename2, LINE, LINE2, LINE3
    character (LEN=4) :: repstr
    character (LEN=16) :: myfmt
    logical :: file_exists
   
    if (errorflag==1) return
   
    allocate (valid(reptot))
    valid = 0
    lines = 0
    if (mod(timeend-timestrt,dtinit).ne.0) timeend=timeend - (mod(timeend-timestrt,dtinit)) + dtinit
    timesteps = (timeend-timestrt)/dtinit + 2
    allocate (output(cols,timesteps), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Allocation error in output array"
      call flush(0)
      errorflag=1
      return
    end if
    
    output = 0.0d0
    tot = 0
   
    allocate(fileline(cols))
   
    do i=1,reptot
      fileunin = 1150 + i
      fileunout = 7150 + i
   
      lines=0
      LINE=""
      
      write(repstr,"(i4.4)") i  
   
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
        write(0,"(a,a,a)") "Error in opening ", trim(filename1), " file"
        errorflag=1
        return
      end if
   
      OPEN(UNIT=fileunout, FILE=trim(filename2),STATUS='UNKNOWN', iostat=ierr)
      if (ierr .ne. 0) then
        write(0,"(a,a,a)") "Error in opening ", trim(filename2), " file"
        errorflag=1
        return
      end if
   
      if (i==1) then
        read(fileunin,"(a141)",iostat=ierr) LINE2
        if (ierr .ne. 0) then
          write(0,"(a,a,a)") "Error in reading from ", trim(filename1), " file to get headers"
          errorflag=1  
          return  
        end if
        rewind (fileunin)
      end if
   
      do while (LINE/="0.00000000E+000")
        read(fileunin,*,iostat=ierr) LINE
        if (ierr .ne. 0) then
          write(0,"(a,a,a)") "Error in reading from ", trim(filename1), " file to skip to data"
          errorflag=1
          return
        end if
      end do
      backspace (fileunin)
     
      do while (ierr==0)
        LINE3=LINE
        lines = lines + 1
        read(fileunin,*,iostat=ierr) LINE
        if ((LINE==LINE3).and.(LINE/="0.00000000E+000")) lines = lines-1      ! Stops any line being re-read
      end do
   
      allocate(input(cols,lines), stat=ierr)
      if (ierr/=0) then
        write(0,"(a)") "Allocation error in input array"
        errorflag=1
        return
      end if
      
      rewind (fileunin)
      do while (LINE/="0.00000000E+000")
        read(fileunin,*,iostat=ierr) LINE
        if (ierr .ne. 0) then
          write(0,"(a,a,a)") "Error in reading from ", trim(filename1), " file to skip to data again"
          errorflag=1
          return
        end if
      end do
      backspace (fileunin)
     
      ierr=0
   
      write(0,*)"total number of lines was ", lines
      do j=1,lines
        read(fileunin,*,iostat=ierr) fileline(1:cols)
        if (ierr/=0) then
          write(0,"(a,a,a,i0)") "Error writing to input array for rep ", repstr, " on line ", j
          errorflag=1
          return
        end if
        do k=1,cols
          input(k,j) = fileline(k)
        end do
      end do
   
      call nevillealg(input,output,cols,lines,errorflag)   ! Interpolate subroutine for single file
      if (errorflag==1) return

      write(myfmt,'(a,i0,a)') '(', cols, '(1x,e16.8e3))'
   
      do k=1,timesteps
        write(fileunout,myfmt) output(:,k)
      end do 
   
      rewind (fileunout)
   
      close (fileunin)
      
      deallocate(input, stat=ierr)
      if (ierr/=0) then
        write(0,"(a)") "Dellocation error in input array"
        errorflag=1
        return
      end if
    end do  
    deallocate(output, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Dellocation error in output array"
      call flush(0)
      errorflag=1
      return
    end if
   
    OPEN(UNIT=501, FILE="normpop.out",STATUS='new', iostat=ierr)
   
    write (501,*), LINE2
    write (501,*), ""
    write (501,*), ""
   
    allocate (output2(cols,timesteps))
   
    output2 = 0.0d0
    fileline = 0.0d0 
   
    do i=1,reptot
   
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
   
    write(myfmt,'(a,i0,a)') '(', cols, '(1x,e16.8e3))'
   
    do k=1,timesteps
      write(501,myfmt) output2(:,k)
    end do

    close (501)
   
    deallocate(valid,output2,fileline, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Deallocation error in valid, output2 and fileline arrays"
      errorflag=1
      return
    end if
   
    return

   end subroutine interpolate

!--------------------------------------------------------------------------------------------------

   subroutine nevillealg (input,output,cols,lines,errorflag)

    implicit none

    real(kind=8), dimension (:,:), intent (in) :: input
    real(kind=8), dimension (:,:), intent (inout) :: output
    integer, intent (inout) :: errorflag
    integer, intent (in) :: cols, lines
    real(kind=8) :: time, timetmp, dataout, errout
    real(kind=8), dimension (:), allocatable :: timein, datain, timeall
    integer :: i, k, n, z, jdown, jup, jold, inpsize

    if (errorflag==1) return

    time = timestrt

    allocate(timeall(lines))

    do i=1,lines
      timeall(i) = input(1,i)
    end do

    n=4

    if (lines.le.n) then
      write(0,"(a)") "File Error! There is not enough data in the files to interpolate!"
      errorflag=1
      return
    end if

    if (input(1,1) == timestrt) then
      output(:,1) = input (:,1)
      time = timestrt
    else
      write(0,"(a)") "Something is wrong. Input array does not start at start time"
      errorflag = 1
      return
    end if
  
    jold = 1
    z=1
    timetmp=timestrt

    do while (time.lt.timeend)

      z=z+1

      if (timeend-dtinit.lt.time) then
        timetmp = timeend
      else
        timetmp = timetmp + dtinit
      end if

      jdown = jold

      call hunt(timeall,lines,timetmp,jdown,errorflag)
      if (errorflag==1) return

      jup = jdown + 1
      jold = jdown

      if (timetmp==timeend) then
        n=3
      end if
      
      inpsize=size(input,2)

      if (mod(n,2)==0) then
        jup = jup + (n/2) - 1
      else
        if (timetmp/=timeend) then
          if ((abs(input(1,min(jup,inpsize))-timetmp)).lt.(abs(input(1,min(jdown,inpsize-1))-timetmp))) then
            jup = jup + (n/2)
          else 
            jup = jup + (n/2) - 1
          end if
        else
          jup = min(jup,inpsize)
        end if
      end if
      
      jup=min(jup,inpsize)
      jdown=jup-n+1

      allocate (timein(n), datain(n))

      do i=1,n
        timein(i) = input(1,jdown+i-1)
      end do
   
      do k=2,cols
        do i=1,n
          datain(i) = input(k,jdown+i-1)
        end do

        call polint (timein, datain, n, timetmp, dataout, errout, errorflag)
        if (errorflag==1) return

        output(k,z) = dataout

      end do
     
      time = timetmp

      output(1,z) = time

      deallocate(timein, datain)

    end do  

   end subroutine nevillealg

!--------------------------------------------------------------------------------------------------

   subroutine polint (xa,ya,n,x,y,dy,errorflag)

    implicit none

    integer, intent (in) :: n
    integer, intent (inout) :: errorflag
    integer :: NMAX, i, m, na
    real(kind=8), intent (inout) :: dy, x, y
    real(kind=8), dimension (:), intent (inout) :: xa, ya
    real(kind=8), dimension (:), allocatable :: c, d
    real(kind=8) :: den, dif, dift, ho, hp, w

    if (errorflag==1) return

    na=1
    NMAX=10

    allocate (c(NMAX), d(NMAX))

    dif = abs(x-xa(1))
    do i=1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
        na=i
        dif=dift
      end if
      c(i)=ya(i)
      d(i)=ya(i)
    end do
    y=ya(na)
    na=na-1
    do m=1,n-1
      do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if (den.eq.0) then  
          write(0,"(3(a,es16.8e3),2(a,i0))") "Error in polint when t=", x, "for xa values", &
                              xa(i)," and", xa(i+m), "when i and m are ",i, "and ", m
          errorflag=1
          return
        end if
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
      end do
      if (2*na.lt.n-m) then
        dy=c(na+1)
      else
        dy=d(na)
        na=na-1
      end if
      y=y+dy
    end do
    return

   end subroutine polint

!-------------------------------------------------------------------------------------------------

   subroutine hunt(timeall,lines,timetmp,jdown,errorflag)

    implicit none

    integer, intent (inout) :: jdown, errorflag
    integer, intent (in) :: lines
    real(kind=8), intent (in) :: timetmp
    real(kind=8), dimension(:), intent (in) :: timeall
    integer :: inc, jup, jm
    logical :: ascnd

    if (errorflag==1) return

    ascnd=timeall(lines).gt.timeall(1)

    if ((timetmp.ge.timeall(jdown)).or.(jdown.gt.lines)) then
      jdown = 0
      jup=lines+1
      goto 3
    end if
    inc=1
    if ((timetmp.ge.timeall(jdown)).eqv.ascnd) then
1      jup=jdown+inc
      if (jup.gt.lines) then
        jup=lines+1
      else if ((timetmp.ge.timeall(jup)).eqv.ascnd) then
        jdown=jup
        inc=inc+inc
        goto 1
      end if
    else
      jup = jdown
2      jdown = jup - inc
      if (jdown.lt.1) then
        jdown = 0
      else if ((timetmp.lt.timeall(jdown)).eqv.ascnd) then
        jup=jdown
        inc=inc+inc
        goto 2
      end if
    end if
3   if ((jup - jdown).eq.1) then
      if (timetmp.eq.timeall(lines)) jdown=lines-1
      if (timetmp.eq.timeall(1)) jdown=1
      return
    end if
    jm=(jup+jdown)/2
    if ((timetmp.ge.timeall(jm)).eqv.ascnd) then
      jdown=jm
    else
      jup=jm
    end if
    goto 3

   end subroutine hunt

!*************************************************************************************************!

END MODULE outputs
