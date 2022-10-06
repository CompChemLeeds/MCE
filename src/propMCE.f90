MODULE propMCE

   use alarrays
   use globvars
   use derivsMCE

!*************************************************************************************************!
!*
!*         Time Propagation (MCE) Module
!*           
!*   Contains subroutines for:
!*
!*      1) Taking a single timestep (top level) and choosing adaptive or static step
!*      2) Controlling the adaptive stepsize, calculating the error and changing the step if necessary
!*      3) Taking a step using the Fourth/Fifth order Runge-Kutta-Cash-Karp system for adaptive steps
!*      4) Taking a step using the Fourth Order Runge-Kutta system for static steps
!*      
!*************************************************************************************************!

contains

!*************************************************************************************************!

  subroutine propstep (bs, dt, dtout, dtfin, time, genflg, timestrt_loc, x, reps)

    implicit none
    type(basisfn), dimension (:), intent (inout) :: bs
    real(kind=8), intent (inout) :: dtout, dtfin
    real(kind=8), intent(inout) :: time, timestrt_loc, dt
    integer, intent(in) :: genflg, x, reps
    
    type(basisfn), dimension (:), allocatable :: dbs_dt1, bserr0, tempbs
    integer::k

    if (errorflag .ne. 0) return
  
    call allocbs(dbs_dt1,size(bs))
    call allocbs(tempbs,size(bs))

    call deriv(bs, dbs_dt1, 1, time, genflg, reps, x)

    ! testing what the different values of the basis correspond to. 
    ! if (x == 1) then 
    !   write(6,*) size(dbs_dt1,2)
    ! endif 


    if (step == "A") then

      call allocbs(bserr0,size(bs))

      do k=1,size(bs)
        bserr0(k)%z(1:ndim) = bs(k)%z(1:ndim) + dbs_dt1(k)%z(1:ndim) + tiny(0.0d0)
        bserr0(k)%d_pes(1:npes) = bs(k)%d_pes(1:npes) + dbs_dt1(k)%d_pes(1:npes) + tiny(0.0d0)
        bserr0(k)%s_pes(1:npes) = bs(k)%s_pes(1:npes) + dbs_dt1(k)%s_pes(1:npes) + tiny(0.0d0)
        bserr0(k)%D_big = bs(k)%D_big + dbs_dt1(k)%D_big + tiny(0.0d0)
      end do

      call rkstpctrl (bs, dbs_dt1, dt, bserr0, dtfin, dtout, tempbs, time, genflg, reps, x)

    else if (step == "S") then

      call rk4 (bs, dbs_dt1, dt, dtfin, dtout, tempbs, time, genflg, x, reps)
  
    end if

    if (errorflag .ne. 0) return
    bs = tempbs

    call deallocbs(dbs_dt1)
    call deallocbs(tempbs)
    if (step=="A") call deallocbs(bserr0)
   

    return

  end subroutine propstep    

!--------------------------------------------------------------------------------------------------

  subroutine rkstpctrl(bs,dbs_dt1,dtprev,bserr0,dtfin,dtnext, tempbs, time, genflg, reps, x)

    implicit none
    type(basisfn), dimension (:), intent (inout) :: tempbs
    type(basisfn), dimension (:), intent (in) :: dbs_dt1, bserr0, bs
    type(basisfn), dimension(:), allocatable :: bserr1
    real(kind=8), intent (inout) :: dtprev, time
    real(kind=8), intent (out) :: dtfin, dtnext
    real(kind=8) :: dt, diff, errmin
    real(kind=8), dimension(:), allocatable :: err1z, err0z   
    integer, intent(in) :: genflg, reps, x
    integer::k, r, adap, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(err1z(size(bs)), err0z(size(bs)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocating the z error arrays in rkstpctrl"
      errorflag=1
    end if

    call allocbs(bserr1,size(bs))

    diff = 1.0d-6

    dt = dtprev

    adap = 0

    do

      call rkck45(bs,dbs_dt1,dt,tempbs,bserr1, time, genflg, reps, x)

      do k=1,size(bs)
        err1z(k) = sum(dble(dconjg(bserr1(k)%z(1:ndim))*(bserr1(k)%z(1:ndim))))
        err0z(k) = sum(dble(dconjg(bserr0(k)%z(1:ndim))*(bserr0(k)%z(1:ndim))))
      end do

      errmin = 1.0d10
      do k=1,size(bs)                   !calculates minimum error ratio over each of the derived values
        do r=1,npes
          errmin = min(abs(bserr0(k)%d_pes(r)/bserr1(k)%d_pes(r)), &
                   abs(bserr0(k)%s_pes(r)/bserr1(k)%s_pes(r)), &
                   errmin)
        end do
        errmin = min(abs(bserr0(k)%D_big/bserr1(k)%D_big),&
                 abs(err0z(k)/err1z(k)),&
                 errmin)
      end do

      errmin = errmin * diff

      if ((errmin.lt.1).and.(adap.eq.0)) then 
        dt = sign(max((0.9*abs(dt)*(errmin**0.25)),dtmin),dt)
        if (dt==dtmin) dtnext=dtmin
        if (dt.lt.dtmin) then
          write(0,"(2(a,es16.8e3),a)") "Error! Somehow, dt is smaller than dtmin, with values of ", &
                                          dt, " and ", dtmin, "respectively"
          errorflag=1 
          return
        end if
        if (time+dt.eq.time) then
          write(0,"(2(a,es16.8e3))") "Error! Underflow in step size at time = ", time, " for dt value of ", dt
          errorflag = 1
          return
        end if
        adap = 1
        cycle
      else if (errmin.ge.1) then
        dtnext = sign(min(0.9*abs(dt)*(errmin**0.2),dtmax),dt)
        exit
      else
        dtnext = sign(min(0.9*abs(dt)*(errmin**0.25),dtmax),dt)
        exit          
      end if

    end do

    dtfin = dt

    call deallocbs(bserr1)
    deallocate(err1z, err0z, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocating the z error arrays in rkstpctrl"
      errorflag=1
    end if

    return

  end subroutine rkstpctrl    

!--------------------------------------------------------------------------------------------------

  subroutine rkck45(bsin,dbs_dt1,dt,tempbs,errbs, time, genflg, reps, x)

    implicit none 
    
    type(basisfn), dimension (:), intent (in) :: bsin, dbs_dt1
    type(basisfn), dimension (:), intent (inout) :: tempbs, errbs
    type(basisfn), dimension (:,:), allocatable::dbs_dt
    real(kind=8),intent(inout)::dt, time
    real(kind=8),dimension(:), allocatable::a,c,d
    real(kind=8),dimension(:,:), allocatable::b
    integer, intent(in) :: genflg, reps, x
    integer::k, l, n, ierr

    if (errorflag .ne. 0) return
    ierr = 0
    allocate (a(6), b(6,6), c(6), d(6), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of Butcher tableau arrays for rkck45"
      errorflag=1
      return
    end if


! Coeffs for Runge-Kutta Cash-Karp method

    a(1)=0.0d0
    a(2)=1.0d0/5.0d0
    a(3)=3.0d0/10.0d0
    a(4)=3.0d0/5.0d0
    a(5)=1.0d0
    a(6)=7.0d0/8.0d0
    b=0.0d0
    b(2,1)=1.0d0/5.0d0
    b(3,1)=3.0d0/40.0d0
    b(3,2)=9.0d0/40.0d0
    b(4,1)=3.0d0/10.0d0
    b(4,2)=-9.0d0/10.0d0
    b(4,3)=6.0d0/5.0d0
    b(5,1)=-11.0d0/54.0d0
    b(5,2)=5.0d0/2.0d0
    b(5,3)=-70.0d0/27.0d0
    b(5,4)=35.0d0/27.0d0
    b(6,1)=1631.0d0/55296.0d0
    b(6,2)=175.0d0/512.0d0
    b(6,3)=575.0d0/13824.0d0
    b(6,4)=44275.0d0/110592.0d0
    b(6,5)=253.0d0/4096.0d0
    c(1)=37.0d0/378.0d0
    c(2)=0.0d0
    c(3)=250.0d0/621.0d0
    c(4)=125.0d0/594.0d0
    c(5)=0.0d0
    c(6)=512.0d0/1771.0d0
    d(1)=c(1)-2825.0d0/27648.0d0
    d(2)=0.0d0
    d(3)=c(3)-18575.0d0/48384.0d0
    d(4)=c(4)-13525.0d0/55296.0d0
    d(5)=c(5)-277.0d0/14336.0d0
    d(6)=c(6)-1.0d0/4.0d0

    allocate (dbs_dt(6,size(bsin)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of RK array"
      errorflag=1
      return
    end if


    do k=1,size(dbs_dt,2)
      do l=1,6
        call allocbf(dbs_dt(l,k))
      end do
    end do

    do k=1,size(dbs_dt,2)
      dbs_dt(1,k)=dbs_dt1(k)
    end do

    do n=2,6

      tempbs = bsin

      do k=1,size(tempbs)
        do l=1,n-1
          tempbs(k)%z(1:ndim)=tempbs(k)%z(1:ndim)+(b(n,l)*dt*dbs_dt(l,k)%z(1:ndim))
          tempbs(k)%d_pes(1:npes)=tempbs(k)%d_pes(1:npes)+(b(n,l)*dt*dbs_dt(l,k)%d_pes(1:npes))
          tempbs(k)%s_pes(1:npes)=tempbs(k)%s_pes(1:npes)+(b(n,l)*dt*dbs_dt(l,k)%s_pes(1:npes))
          tempbs(k)%a_pes(1:npes)=tempbs(k)%d_pes(1:npes)*exp(i*tempbs(k)%s_pes(1:npes))
          tempbs(k)%D_big=tempbs(k)%D_big+(b(n,l)*dt*dbs_dt(l,k)%D_big)
        end do
      end do

      call deriv(tempbs, dbs_dt(n,:), n, time, genflg, reps, x)

    end do

    tempbs = bsin

    do k=1,size(tempbs)
      do l=1,6
        tempbs(k)%z(1:ndim)=tempbs(k)%z(1:ndim)+(c(l)*dt*dbs_dt(l,k)%z(1:ndim))
        tempbs(k)%d_pes(1:npes)=tempbs(k)%d_pes(1:npes)+(c(l)*dt*dbs_dt(l,k)%d_pes(1:npes))
        tempbs(k)%s_pes(1:npes)=tempbs(k)%s_pes(1:npes)+(c(l)*dt*dbs_dt(l,k)%s_pes(1:npes))
        tempbs(k)%a_pes(1:npes)=tempbs(k)%d_pes(1:npes)*exp(i*tempbs(k)%s_pes(1:npes))
        tempbs(k)%D_big=tempbs(k)%D_big+(c(l)*dt*dbs_dt(l,k)%D_big)
      end do
    end do

    do k=1,size(tempbs)
      do l=1,6
        errbs(k)%z(1:ndim)=errbs(k)%z(1:ndim)+(d(l)*dbs_dt(l,k)%z(1:ndim))
        errbs(k)%d_pes(1:npes)=errbs(k)%d_pes(1:npes)+(d(l)*dbs_dt(l,k)%d_pes(1:npes))
        errbs(k)%s_pes(1:npes)=errbs(k)%s_pes(1:npes)+(d(l)*dbs_dt(l,k)%s_pes(1:npes))
        errbs(k)%D_big=errbs(k)%D_big+(d(l)*dbs_dt(l,k)%D_big)
      end do
    end do

    do k=1,size(dbs_dt,2)
      do l=1,size(dbs_dt,1)
        call deallocbf(dbs_dt(l,k))
      end do
    end do    
    deallocate (a,b,c,d,dbs_dt, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of arrays or matrices for rkck45"
      errorflag=1
      return
    end if

    return

  end subroutine rkck45

!--------------------------------------------------------------------------------------------------

  subroutine rk4 (bsin, dbs_dt1, dt, dtfin, dtout, tempbs, time, genflg, x, reps)

    implicit none 
    
    type(basisfn), dimension (:), intent (in) :: bsin, dbs_dt1
    type(basisfn), dimension (:), intent (inout) :: tempbs
    type(basisfn), dimension (:,:), allocatable :: dbs_dt
    real(kind=8), dimension(ndim) :: dummy_arr
    real(kind=8),intent(inout)::dtfin, dtout
    real(kind=8), intent(inout) :: time, dt
    real(kind=8),dimension(:), allocatable::a,b,c
    integer, intent(in) :: genflg, x, reps
    integer::k, l, n, ierr

    if (errorflag .ne. 0) return
    ierr=0
    allocate(a(4), b(4), c(4), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of rk4 Butcher tableau parameters"
      errorflag=1
      return
    end if 
    
    dummy_arr=0.0d0

! Coeffs for Runge-Kutta 4 method

    a(1)=0.0d0
    a(2)=1.0d0/2.0d0
    a(3)=1.0d0/2.0d0
    a(4)=1.0d0
    b(1)=0.0d0
    b(2)=1.0d0/2.0d0
    b(3)=1.0d0/2.0d0
    b(4)=1.0d0
    c(1)=1.0d0/6.0d0
    c(2)=1.0d0/3.0d0
    c(3)=1.0d0/3.0d0
    c(4)=1.0d0/6.0d0

    allocate (dbs_dt(4,size(bsin)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in basis set allocation in rk4"
      errorflag=1
      return
    end if

    do k=1,size(dbs_dt,2)
      do l=1,size(dbs_dt,1)
        call allocbf(dbs_dt(l,k))
      end do
    end do
    do k=1,size(dbs_dt,2)
      dbs_dt(1,k)=dbs_dt1(k)
    end do
    
    do n=2,4

      do k=1,size(tempbs)
        tempbs(k)%z(1:ndim)=bsin(k)%z(1:ndim)+(b(n)*dt*dbs_dt(n-1,k)%z(1:ndim))
        tempbs(k)%d_pes(1:npes)=bsin(k)%d_pes(1:npes)+(b(n)*dt*dbs_dt(n-1,k)%d_pes(1:npes))
        tempbs(k)%s_pes(1:npes)=bsin(k)%s_pes(1:npes)+(b(n)*dt*dbs_dt(n-1,k)%s_pes(1:npes))
        tempbs(k)%D_big=bsin(k)%D_big+(b(n)*dt*dbs_dt(n-1,k)%D_big)
        tempbs(k)%a_pes(1:npes)=tempbs(k)%d_pes(1:npes)*exp(i*tempbs(k)%s_pes(1:npes))
      end do

      call deriv(tempbs, dbs_dt(n,:), n, time, genflg, reps, x)

    end do

    tempbs=bsin

    do k=1,size(tempbs)
      do l=1,4
        tempbs(k)%z(1:ndim)=tempbs(k)%z(1:ndim)+(c(l)*dt*dbs_dt(l,k)%z(1:ndim))
        tempbs(k)%d_pes(1:npes)=tempbs(k)%d_pes(1:npes)+(c(l)*dt*dbs_dt(l,k)%d_pes(1:npes))
        tempbs(k)%s_pes(1:npes)=tempbs(k)%s_pes(1:npes)+(c(l)*dt*dbs_dt(l,k)%s_pes(1:npes))
        tempbs(k)%a_pes(1:npes)=tempbs(k)%d_pes(1:npes)*exp(i*tempbs(k)%s_pes(1:npes))
        tempbs(k)%D_big=tempbs(k)%D_big+(c(l)*dt*dbs_dt(l,k)%D_big)
      end do
    end do

    dtfin = dt
    dtout = dt

    do k=1,size(dbs_dt,2)
      do l=1,size(dbs_dt,1)
        call deallocbf(dbs_dt(l,k))
      end do
    end do
    deallocate (dbs_dt, a,b,c, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of local arrays in rk4"
      errorflag=1
      return
    end if

    return

  end subroutine rk4

!*************************************************************************************************!

end module propMCE

