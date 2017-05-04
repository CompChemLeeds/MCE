MODULE bsetgen
  
  use globvars
  use Ham
  use alarrays
  use outputs
  use Chks
  use redirect
  use propMCE

!***********************************************************************************!
!*
!*         Basis Set generation Module
!*           
!*   Contains subroutines for:
!*
!*      1) Redirecting to the appropriate basis set generation subroutine
!*      2) Generating the basis set coherent states with a swarm structure
!*      3) Generating the basis set coherent states with a train structure
!*      4) Generating the basis set coherent states with a swarm/train structure
!*      5) Generating the basis set coherent states with a grid structure
!*      6) Generating the basis set coherent states with a grid/swarm structure
!*      7) Generating the basis set multiconfigurational D prefactor
!*      8) Generating the basis function single configurational d prefactors for 
!*              electronic states
!*      
!***********************************************************************************!

contains

!***********************************************************************************!
!              Populate/Generate array Subroutines
!***********************************************************************************!

  subroutine genbasis(bs, mup, muq, alcmprss, gridsp, t, initgrid, reps, map_bfs) 

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    complex(kind=8), dimension(:,:), intent(inout)::initgrid
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(inout) :: alcmprss, gridsp, t
    integer, intent(inout), dimension(:,:) :: map_bfs
    integer, intent(in)::reps

    if (errorflag .ne. 0) return

    select case (basis)
      case ("SWARM")
        call gen_swarm(bs,mup,muq,alcmprss,t)
      case ("TRAIN")
        call gen_train(bs,mup,muq,t,reps,map_bfs)
      case ("SWTRN")
        call gen_swtrn(bs,mup,muq,alcmprss,t,reps,map_bfs)
      case ("GRID")
        call gen_grid(bs,mup,muq,initgrid,gridsp)
      case ("GRSWM")
        call gen_grswm(bs,mup,muq,alcmprss,initgrid,gridsp,t)
      case default
        write(0,"(a)") "Error! Initial basis calculation method not recognised!"
        write(0,"(a)") "This should have been caught at the read stage!"
        errorflag = 1
        return
    end select

    return

  end subroutine genbasis

!------------------------------------------------------------------------------------

  subroutine gen_swarm(bs,mup,muq,alcmprss,t)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: alcmprss, t
    type (basisfn) :: bf
    integer::m, k, n, h, ierr, redo

    if (errorflag .ne. 0) return
    ierr = 0
    n=0

    call allocbf(bf)
    
    
    if (sys.eq."HP") then
!    if (.false.) then
    
    h=0.5*size(bs)
    
    do k=1,h
      bf=bs(k)
      do
        do m=1,ndim
!          bf%z(m)=gauss_random(alcmprss,muq(m),mup(m))
          bf%z(m)=cmplx((ZBQLNOR(muq(m)+5.0d0,sigq*sqrt(alcmprss))) &
         ,((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*sqrt(alcmprss)))),kind=8)
        end do
        call enchk(bf,t,n,redo,k)
        if (redo==1) cycle
        if (redo==0) exit
      end do
      bs(k)=bf    
    end do  
    do k=h+1, size(bs)
      bf=bs(k)
      do
        do m=1,ndim
!          bf%z(m)=gauss_random(alcmprss,muq(m),mup(m))
          bf%z(m)=cmplx((ZBQLNOR(muq(m)-5.0d0,sigq*sqrt(alcmprss))) &
         ,((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*sqrt(alcmprss)))),kind=8)
        end do
        call enchk(bf,t,n,redo,k)
        if (redo==1) cycle
        if (redo==0) exit
      end do
      bs(k)=bf    
    end do  
    
    else
    
    do k=1,size(bs)
      bf=bs(k)
      do
        do m=1,ndim
!          bf%z(m)=gauss_random(alcmprss,muq(m),mup(m))
          bf%z(m)=cmplx((ZBQLNOR(muq(m),sigq*sqrt(alcmprss))) &
         ,((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*sqrt(alcmprss)))),kind=8)
        end do
        call enchk(bf,t,n,redo,k)
        if (redo==1) cycle
        if (redo==0) exit
      end do
      bs(k)=bf
    end do
    
    end if

    call deallocbf(bf)

  end subroutine gen_swarm

!------------------------------------------------------------------------------------

  subroutine gen_train(bs,mup,muq,t,reps,map_bfs)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(inout) :: t
    integer, intent(inout), dimension(:,:) :: map_bfs
    integer, intent(in)::reps
    type (basisfn), dimension(:), allocatable :: bf
    real(kind=8) :: dt, dtnext, dtdone, timeold
    integer::m, n, ierr, stepback, x, genflg

    if (errorflag .ne. 0) return
    ierr=0
    n=0
    genflg = 1

    timeold = t

    call allocbs(bf,1)

    do m=1,ndim
      bf(1)%z(m) = cmplx(muq(m),mup(m),kind=8)
    end do

    bf(1)%D_big = (1.0d0,0.0d0)
    bf(1)%d_pes(in_pes) = (1.0d0,0.0d0)
    bf(1)%a_pes(in_pes) = (1.0d0,0.0d0)
    
    if (mod(in_nbf,2)==1) stepback = ((in_nbf-1)/2)*trainsp
    if (mod(in_nbf,2)==0) stepback = ((in_nbf/2)*trainsp)-(trainsp/2)
    dt = -1.0d0*dtinit

    do x=1,stepback
      call propstep(bf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps,map_bfs)
      t = t + dt
    end do

    write(6,"(a,es16.8e3)") "t at maximum stepback is ", t

    dt = dtinit

    do x=1,(((in_nbf-1)*trainsp)+1)
      if (mod(x-1,trainsp)==0) then
        bs(((x-1)/trainsp)+1)%z = bf(1)%z
        bs(((x-1)/trainsp)+1)%a_pes = bf(1)%a_pes
        bs(((x-1)/trainsp)+1)%d_pes = bf(1)%d_pes
        bs(((x-1)/trainsp)+1)%s_pes = bf(1)%s_pes
      end if
      call propstep(bf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps,map_bfs)
      t = t + dt
    end do

    write(6,"(a,es16.8e3)") "t at maximum stepforward is ", t

    call deallocbs(bf)

    t = timeold

    return       

  end subroutine gen_train

!------------------------------------------------------------------------------------

  subroutine gen_swtrn(bs,mup,muq,alcmprss,t,reps,map_bfs)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(inout) :: alcmprss, t
    integer, intent(inout), dimension(:,:) :: map_bfs
    integer, intent(in)::reps
    type(basisfn), dimension(:), allocatable :: swrmbf, tmpbf
    type(basisfn) :: bf
    real(kind=8) :: dt, dtdone, dtnext, timeold
    integer::m, k, j, x, n, ierr, trtype, swrmsize, restart, kcut, kstrt, steps, redo
    integer::def_stp2, stepback, genflg

    if (errorflag .ne. 0) return
    ierr=0
    n=0
    restart = 0
    trtype = 1 
    def_stp2 = def_stp
    genflg = 1

    if (mod(def_stp2,2)==0) def_stp2 = def_stp2 + 1

    if (mod(in_nbf,def_stp2)==0) then
      steps = def_stp2
    else if (mod(in_nbf,def_stp2 - 1)==0) then
      steps = def_stp2 - 1
    else
      write(0,'(a,i0,a,i0)'), "In_nbf is not an integer multiple of ", def_stp2, &
                 " or ", def_stp2 - 1, ", which are the default number of steps."
      write(0,'(a,a)'), "To have a different number of steps in a swarm train or ",&
                    "train swarm, modify value of def_stp in gen_swtrn sub" 
      errorflag = 1
      return
    end if

    swrmsize = in_nbf / steps

    write(6,"(a,i0,a,i0,a)")"Making ",swrmsize," trains each ",steps," carriages long."
    
    if ((mod(swrmsize,2)==1).and.(mod(steps,2)==1)) uplimnorm = 1.0000001d0

    timeold = t

    call allocbs(swrmbf,swrmsize)
    call allocbf(bf)

    if (mod(swrmsize,2)==0) then
      kcut = 1
    else
      do m=1,ndim
        swrmbf(1)%z(m) = cmplx(muq(m),mup(m),kind=8)
      end do
      swrmbf(1)%d_pes(in_pes) = (1.0d0,0.0d0)
      swrmbf(1)%a_pes(in_pes) = (1.0d0,0.0d0)
      kcut = 2
    end if

    do k=kcut,size(swrmbf)
      do
        do m=1,ndim
          bf%z(m)=gauss_random(alcmprss,muq(m),mup(m))
!          bf%z(m)=cmplx((ZBQLNOR(muq(m),sigq*sqrt(alcmprss))) &
!         ,((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*sqrt(alcmprss)))),kind=8)
        end do
        call enchk(bf,t,n,redo,k)
        if (redo==1) cycle
        if (redo==0) exit
      end do
      bf%d_pes(in_pes) = (1.0d0,0.0d0)
      bf%a_pes(in_pes) = (1.0d0,0.0d0)
      swrmbf(k)=bf
    end do

    if (trtype == 1) then

      do k=1,size(swrmbf)
        swrmbf%D_big = (1.0d0,0.0d0)
      end do
    
      if (mod(steps,2)==1) stepback = ((steps-1)/2)*trainsp
      if (mod(steps,2)==0) stepback = ((steps/2)*trainsp)-(trainsp/2)
      dt = -1.0d0*dtinit

      do x=1,stepback
        call propstep(swrmbf,dt,dtnext,dtdone,t,genflg,timestrt,x-stepback,reps,map_bfs)
        t = t + dt
      end do

      dt = dtinit

      kcut = steps*trainsp     

      do x=1,kcut
        if (mod(x-1,trainsp)==0) then
          do j=1,swrmsize
            bs((((x-1)*swrmsize)/trainsp)+j)%z = swrmbf(j)%z
            bs((((x-1)*swrmsize)/trainsp)+j)%a_pes = swrmbf(j)%a_pes
            bs((((x-1)*swrmsize)/trainsp)+j)%d_pes = swrmbf(j)%d_pes
            bs((((x-1)*swrmsize)/trainsp)+j)%s_pes = swrmbf(j)%s_pes
          end do
        end if
        call propstep(swrmbf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps,map_bfs)
        t = t + dt
      end do

    else if (trtype==2) then

      call allocbs(tmpbf,1)

      do j=1,swrmsize

        tmpbf(1) = swrmbf(j)

        tmpbf%D_big = (1.0d0,0.0d0)
    
        if (mod(steps,2)==1) stepback = ((steps-1)/2)*trainsp
        if (mod(steps,2)==0) stepback = ((steps/2)*trainsp)-(trainsp/2)
        dt = -1.0d0*dtinit

        do x=1,stepback
          call propstep(tmpbf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps,map_bfs)
          t = t + dt
        end do

        dt = dtinit

        kstrt = (steps*(j-1)*trainsp)+1
        kcut = steps * j * trainsp   

        do x=kstrt,kcut
          if (mod(x-1,trainsp)==0) then
            bs(((x-1)/trainsp)+1)%z = tmpbf(1)%z
            bs(((x-1)/trainsp)+1)%a_pes = tmpbf(1)%a_pes
            bs(((x-1)/trainsp)+1)%d_pes = tmpbf(1)%d_pes
            bs(((x-1)/trainsp)+1)%s_pes = tmpbf(1)%s_pes
          end if
          call propstep(tmpbf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps,map_bfs)
          t = t + dt
        end do

      end do

      call deallocbs(tmpbf)

    else

      write(0,'(a,a,i0)'), "Error! Train swarm/swarm train identifier is wrong. ",&
                 "Should be 1 or 2 but got ", trtype
      errorflag = 1
      return

    end if       

    call deallocbs(swrmbf)
    call deallocbf(bf)

    t = timeold

  end subroutine gen_swtrn

!------------------------------------------------------------------------------------
  subroutine gen_grid(bs,mup,muq,initgrid,gridsp)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    complex(kind=8), dimension(:,:), intent(inout)::initgrid
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: gridsp
    type (basisfn) :: bf
    real(kind=8), dimension (:), allocatable ::qstrt, pstrt
    real(kind=8) :: realz, imagz
    integer, dimension (:,:), allocatable :: coordinates
    integer, dimension (:), allocatable :: qsize, psize, columnsize, columnrep
    integer::m, k, l, n, x, y, ierr

    if (errorflag .ne. 0) return
    ierr=0
    n=0

    k=0
    allocate (qstrt(ndim), stat=ierr)
    if (ierr==0) allocate (pstrt(ndim), stat=ierr)
    if (ierr==0) allocate (qsize(ndim), stat=ierr)
    if (ierr==0) allocate (psize(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation starting p and q values for grid"
      errorflag=1
      return
    end if

    if (ndim==3) then
      qsize(1) = max(qsizez,1)
      qsize(2) = max(qsizex,1)
      qsize(3) = max(qsizey,1) 
      psize(1) = max(psizez,1)
      psize(2) = max(psizex,1)
      psize(3) = max(psizey,1) 
    else if (ndim==1) then
      qsize(1) = qsizez
      psize(1) = psizez
    else
      write(0,"(a)") "Using a grid of neither 3 nor 1 dimensions! How did that happen?"
      errorflag=1
      return
    end if

    allocate (coordinates(in_nbf,ndim*2))
    allocate (columnsize(ndim*2))
    allocate (columnrep(ndim*2))

    n=0
    do m=1,ndim*2
      if (mod(m,2)==1) then
        n=n+1
        columnsize(m)=qsize(n)
      else
        columnsize(m)=psize(n)
      end if
      columnrep(m) = 1
      if (m/=1) then
        do l=1,m-1
          columnrep(m) = columnrep(m)*columnsize(l)
        end do
      end if
    end do
    do m=1,ndim*2
      columnrep(m) = 1
      if (m/=1) then
        do l=1,m-1
          columnrep(m) = columnrep(m)*columnsize(l)
        end do
      end if
    end do

    do m=1,ndim*2
      x=in_nbf/(columnrep(m)*columnsize(m))
      y=in_nbf/columnrep(m)
      do n=1,in_nbf
        if (n==1) then
          coordinates(n,m) = 1
        else if (mod(n-1,x)==0) then 
          if (mod(n-1,y)==0) then
            coordinates(n,m) = 1
          else
            coordinates(n,m) = coordinates(n-1,m) + 1
          end if
        else
          coordinates(n,m) = coordinates(n-1,m)
        end if
      end do
    end do

    do m=1,ndim
      qstrt(m) = muq(m)-((dble(qsize(m))-1.0d0)*sigq*gridsp/2.0d0)
      pstrt(m) = mup(m)-((dble(psize(m))-1.0d0)*sigp*gridsp/2.0d0)
    end do

    do k=1,in_nbf
      bf=bs(k)
      do m=1,ndim
        realz = qstrt(m)+(gridsp*((coordinates(k,(2*m)-1)-1)*sigq))
        imagz = pstrt(m)+(gridsp*((coordinates(k,(2*m)  )-1)*sigp))
        bf%z(m)=cmplx(realz, imagz,kind=8) 
      end do
      do m=1,ndim
        initgrid(k,m)=bf%z(m)  
      end do
      bs(k)=bf
    end do
    
    if ((mod(qsizez,2)==0).and.(mod(in_nbf,2)==1)) then
      do m=1,ndim
        bf%z(m)=cmplx(muq(m),mup(m),kind=8)
      end do
      do m=1,ndim
        initgrid(in_nbf,m)=bf%z(m)
      end do
      bs(in_nbf)=bf
    end if

    if (size(initgrid,1).ne.in_nbf) then
      write(0,"(a,a)") "Error in generating the initial grid overlap matrix. ",&
                "Sizes are not equal to in_nbf"
      write(0,"(a,i0)") "Error is in size of initgrid, which is ", size(initgrid)
      write(0,"(a,i0)") "in_nbf is ", in_nbf          
      errorflag = 1
      return
    end if

    deallocate (qstrt, stat=ierr)
    if (ierr==0) deallocate (pstrt, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation starting p and q values for grid"
      errorflag=1
      return
    end if
    call deallocbf(bf)

  end subroutine gen_grid

!------------------------------------------------------------------------------------

  subroutine gen_grswm(bs,mup,muq,alcmprss,initgrid,gridsp,t)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    complex(kind=8), dimension(:,:), intent(inout)::initgrid
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: alcmprss, gridsp, t
    type (basisfn) :: bf
    real(kind=8):: qstrt, pstrt
    integer::m, k, p, q, n, ierr, redo

    if (errorflag .ne. 0) return
    ierr=0
    n=0

    k=0

    qstrt = muq(1)-((dble(qsizez)-1.0d0)*sigq*gridsp/2.0d0)
    pstrt = mup(1)-((dble(psizez)-1.0d0)*sigp*gridsp/2.0d0)
    do p=1,psizez
      do q=1,qsizez
        k=q+((p-1)*qsizez)
        bf=bs(k)
        do
          bf%z(1)=cmplx((qstrt+(gridsp*(q-1)*sigq)),&
                        (pstrt+(gridsp*(p-1)*sigp)),kind=8)
          do m=2,ndim
            bf%z(m)=gauss_random(1.0d0/alcmprss,muq(m),mup(m))
!            bf%z(m)=cmplx((ZBQLNOR(muq(m),sigq*sqrt(alcmprss))),&
!            ((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*sqrt(alcmprss)))),kind=8)
          end do
          call enchk(bf,t,n,redo,k)
          if (redo==1) cycle
          if (redo==0) exit
        end do
        do m=1,ndim
          initgrid(k,m)=bf%z(m)  
        end do
        bs(k)=bf
      end do
    end do
    if ((mod(qsizez,2)==0).and.(mod(in_nbf,2)==1)) then
      do m=1,ndim
        bf%z(m)=cmplx(muq(m),mup(m),kind=8)
      end do
      do m=1,ndim
        initgrid(in_nbf,m)=bf%z(m)
      end do
      bs(in_nbf)=bf
    end if
    if (size(initgrid,1).ne.in_nbf) then
      write(0,"(a,a)") "Error in generating the initial grid overlap matrix. ",&
                "Sizes are not equal to in_nbf"
      write(0,"(a,i0)") "Error is in size of initgrid, which is ", size(initgrid)
      write(0,"(a,i0)") "in_nbf is ", in_nbf          
      errorflag = 1
      return
    end if

    call deallocbf(bf)

  end subroutine gen_grswm

!------------------------------------------------------------------------------------

  subroutine genD_big(bs, mup, muq, restart)   !   Level 1 Subroutine

    implicit none
    type(basisfn),dimension(:),intent(inout)::bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    complex(kind=8),dimension(:,:),allocatable::ovrlp_mat
    complex(kind=8),dimension(:), allocatable::zinit, zpq
    complex(kind=8),dimension(:), allocatable:: C_k, D
    integer, intent(inout) :: restart
    integer::k, j, ierr

    if (errorflag .ne. 0) return
    ierr = 0

    allocate(D(size(bs)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of D_k array in genD_big"
      errorflag=1
      return
    end if 

    if ((basis.ne."TRAIN").and.(basis.ne."SWTRN")) then
      do j=1,size(bs)
        bs(j)%a_pes(in_pes) = (1.0d0,0.0d0)
        bs(j)%d_pes(in_pes) = (1.0d0,0.0d0)
      end do
    end if

    if (((basis=="GRID").or.(basis=="GRSWM")).and.(mod(in_nbf,2)==1)) then

      if (mod(psizez,2)==1) then
        if (basis=="GRSWM") then 
          k = ((qsizez-1)/2)+(((psizez-1)/2)*qsizez)
        else
          k = (in_nbf+1)/2
        end if
      else
        k = in_nbf
      end if
      do j=1,size(bs)
        if (j.ne.k) then
          D(k)=(0.0d0,0.0d0)
        else
          D(k)=(1.0d0,0.0d0)
        end if
      end do

    else

      allocate(C_k(size(bs)), stat = ierr)
      if (ierr==0) allocate(zinit(ndim), stat = ierr)
      if (ierr==0) allocate(zpq(ndim), stat = ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error in allocation of C_k or z arrays in genD_big"
        errorflag=1
        return
      end if

      zinit(1:ndim) = cmplx(muq(1:ndim),mup(1:ndim),kind=8)

      do k=1,size(bs)
        zpq(1:ndim) = bs(k)%z(1:ndim)
        C_k(k) = ovrlpij(zpq, zinit) * dconjg(bs(k)%d_pes(in_pes)) * &
                    cdexp (-1.0d0*i*bs(k)%s_pes(in_pes))
      end do
   
      deallocate(zpq, stat = ierr)
      if (ierr==0) deallocate(zinit, stat = ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error in z array deallocation in genD_big"
        errorflag=1
        return
      end if

      allocate(ovrlp_mat(size(bs),size(bs)), stat = ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error in overlap matrix allocation in genD_big"
        errorflag=1
        return
      end if
      ovrlp_mat=ovrlpphimat(bs)

      if (matfun.eq.'zgesv') then
        call lineq(ovrlp_mat, C_k, D)
      else if (matfun.eq.'zheev') then
        call matinv(ovrlp_mat, C_k, D, restart)
      else
        write(0,"(a)") "Error! Matrix function not recognised! Value is ", matfun
        errorflag = 1
        return
      end if

      deallocate(C_k, stat = ierr)
      if (ierr==0) deallocate(ovrlp_mat, stat = ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error in deallocation of local arrays in genD_big"
        errorflag=1
        return
      end if

    end if

    do k=1,size(bs)
      bs(k)%D_big = D(k)
    end do

    if (method == "MCEv1") then
      do j=1,size(bs)
        bs(j)%a_pes(in_pes) = bs(j)%D_big
        bs(j)%d_pes(in_pes) = bs(j)%D_big
        bs(j)%D_big = (1.0d0,0.0d0)
      end do
    end if

    return

  end subroutine genD_big

!------------------------------------------------------------------------------------

  subroutine gend_small(bs)   !   Level 1 Subroutine

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    integer :: j

    if (errorflag .ne. 0) return

    if (method == "MCEv1") then
      do j=1,size(bs)
        bs(j)%a_pes(in_pes) = bs(j)%D_big
        bs(j)%d_pes(in_pes) = bs(j)%D_big
        bs(j)%D_big = (1.0d0,0.0d0)
      end do
    else if ((method == "MCEv2").or.(method == "CCS")) then
      do j=1,size(bs)
        bs(j)%a_pes(in_pes) = (1.0d0,0.0d0)
        bs(j)%d_pes(in_pes) = (1.0d0,0.0d0)
      end do
    else
      write(0,"(a)") "Error! The method is wrong"
      write(0,"(a)") "I don't even know how you got this far"
      errorflag = 1
      return
    end if

    return

  end subroutine gend_small 

!*************************************************************************************************!

  function gauss_random (width, muq, mup)
  
    implicit none
    complex(kind=8) :: gauss_random
    real(kind=8), intent(in) :: width, mup, muq
    integer :: iset
    real(kind=8) :: fac,rsq,v1,v2,ran_1,xz,yz

    iset = 0
    
    do while (iset .lt. 2)      
      CALL RANDOM_NUMBER(ran_1) 
      v1=2.0d0*ran_1-1.0d0
      CALL RANDOM_NUMBER(ran_1)     
      v2=2.0d0*ran_1-1.0d0
      rsq=v1**2.0d0+v2**2.0d0
      if (rsq.ge.1.0d0.or.rsq.eq.0.0d0) cycle
      fac=sqrt(-2.0d0*log(rsq)/rsq)
      if (iset == 0) then
        xz=v2*fac*sigq
      else if (iset == 1) then
        yz=v2*fac*sigp
      end if    
      iset=iset + 1
    end do
    
    xz = muq + (xz * width)
    yz = mup + (yz * width)
    
    gauss_random = cmplx(xz,yz,kind=8)
    
    return
    
  end function gauss_random
end module bsetgen
