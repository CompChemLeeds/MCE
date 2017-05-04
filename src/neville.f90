module neville

contains

subroutine nevillealg (input,output,dt,timeend,timestrt,cols,lines,errorflag)

  implicit none

  real(kind=8), dimension (:,:), intent (in) :: input
  real(kind=8), dimension (:,:), intent (inout) :: output
  real(kind=8), intent (in) :: dt, timestrt, timeend
  integer, intent (inout) :: errorflag
  integer, intent (in) :: cols, lines
  real(kind=8) :: time, timetmp, dataout, errout
  real(kind=8), dimension (:), allocatable :: timein, datain, timeall
  integer :: i, k, n, z, jdown, jup, jold

  if (errorflag==1) return

  time = timestrt

  allocate(timeall(lines))

  do i=1,lines
     timeall(i) = input(1,i)
  end do

  n=4

  if (lines.le.n) then
     print *, "File Error! There is not enough data in the files to interpolate!"
     errorflag=1
     return
  end if

  if (input(1,1) == timestrt) then
     output(:,1) = input (:,1)
     time = timestrt
  else
     print *, "Something is wrong. Input array does not start at start time"
     errorflag = 1
     return
  end if
 
  if (input(1,2) == timestrt + dt) then
     output(:,2) = input (:,2)
     time = time + dt
  else
     print *, "I expected the first timestep to be the same as the initial timestep"
     errorflag = 1
     return
  end if

  jold = 2
  z=2
  timetmp=time  

  do while (time.lt.timeend)

     z=z+1

     if (timeend-dt.lt.time) then
        timetmp = timeend
     else
        timetmp = timetmp + dt
     end if

     jdown = jold

     call hunt(timeall,lines,timetmp,jdown,errorflag)
     if (errorflag==1) return

     jup = jdown + 1

     if (timetmp==timeend) then
        n=3
     end if

     jold = jdown

     if (mod(n,2)==0) then
        jdown = jdown - (n/2) + 1
        jup = jup + (n/2) - 1
     else
        if (((abs(input(1,jup)-timetmp)).lt.(abs(input(1,jdown)-timetmp))).and.(timetmp/=timeend)) then
           jup = jup + (n/2)
           jdown = jdown - (n/2) + 1
        else
           jup = jup + (n/2) - 1
           jdown = jdown - (n/2)
        end if
     end if

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
           print *, "Error in polint"
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
1    jup=jdown+inc
     if (jup.gt.lines) then
        jup=lines+1
     else if ((timetmp.ge.timeall(jup)).eqv.ascnd) then
        jdown=jup
        inc=inc+inc
        goto 1
     end if
  else
     jup = jdown
2    jdown = jup - inc
     if (jdown.lt.1) then
        jdown = 0
     else if ((timetmp.lt.timeall(jdown)).eqv.ascnd) then
        jup=jdown
        inc=inc+inc
        goto 2
     end if
  end if
3 if ((jup - jdown).eq.1) then
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

end module
