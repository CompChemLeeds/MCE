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

  subroutine genbasis(bs, mup, muq, alcmprss, t, reps) 

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(inout) :: alcmprss, t
    integer, intent(in)::reps

    if (errorflag .ne. 0) return

    select case (basis)
      case ("SWARM")
        call gen_swarm(bs,mup,muq,alcmprss,t,reps)
      case ("SWTRN")
        call gen_swtrn(bs,mup,muq,alcmprss,t,reps)
      case default
        write(0,"(a)") "Error! Initial basis calculation method not recognised!"
        write(0,"(a)") "This should have been caught at the read stage!"
        errorflag = 1
        return
    end select

    return

  end subroutine genbasis

!------------------------------------------------------------------------------------

  subroutine gen_swarm(bs,mup,muq,alcmprss,t,reps)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: alcmprss, t
    integer, intent(in) :: reps
    
    type (basisfn) :: bf
    integer::m, k, n, h, ierr, redo

    if (errorflag .ne. 0) return
    ierr = 0
    n=0

    call allocbf(bf)
      
    do k=1,size(bs)
      bf=bs(k)
      do
        do m=1,ndim
          if (randfunc.eq."GAUS") then
            bf%z(m)=gauss_random(alcmprss,muq(m),mup(m))
          else
            bf%z(m)=cmplx((ZBQLNOR(muq(m),sigq*alcmprss)) &
           ,((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*alcmprss))),kind=8)
          end if
        end do
        call enchk(bf,t,n,redo,k,reps)
        if (redo==1) cycle
        if (redo==0) exit
      end do
      bs(k)=bf
    end do

    call deallocbf(bf)

  end subroutine gen_swarm

!------------------------------------------------------------------------------------

  subroutine gen_swtrn(bs,mup,muq,alcmprss,t,reps)

    implicit none
    type(basisfn), dimension(:), intent(inout) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(inout) :: alcmprss, t
    integer, intent(in)::reps
    type(basisfn), dimension(:), allocatable :: swrmbf, tmpbf
    type(basisfn) :: bf
    real(kind=8) :: dt, dtdone, dtnext, timeold, q, p, sumamps
    integer::m, k, j, x, n, r, ierr, trtype, swrmsize, restart, kcut, kstrt, steps, redo
    integer::def_stp2, stepback, genflg

    if (errorflag .ne. 0) return
    ierr=0
    n=0
    restart = 0
    trtype = 1 
    def_stp2 = train_len
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
                    "train swarm, modify value of def_stp in the input file" 
      errorflag = 1
      return
    end if

    swrmsize = swtrn_swrm    

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
          if (randfunc.eq."GAUS") then
            bf%z(m)=gauss_random(alcmprss,muq(m),mup(m))
          else
            bf%z(m)=cmplx((ZBQLNOR(muq(m),sigq*alcmprss)) &
            ,((1.0d0/hbar)*(ZBQLNOR(mup(m),sigp*alcmprss))),kind=8)
          end if
        end do
        call enchk(bf,t,n,redo,k,reps)
        if (redo==1) cycle
        if (redo==0) exit
      end do
      if (qss==1) then
        if (npes==2) then
          call random_number(q)
          q = q * 2.0 * pirl
          call random_number(p)
          p = p * 2.0 * pirl
          do r=1,npes
            if (r==in_pes) then
              bf%a_pes(r)=cmplx(dcos(q),0.0d0)
            else
              bf%a_pes(r)=cmplx(dsin(q)*dcos(p),dsin(q)*dsin(p))
            end if
          end do
        else
          sumamps = 0.0d0
          do r=1,npes
            call random_number(q)
            call random_number(p)
            q = q * 2.0 - 1.0
            p = p * 2.0 - 1.0
            bf%a_pes(r) = cmplx(q,p)
            sumamps = sumamps + dsqrt(dble(cmplx(q,p) * cmplx(q,-1.0*p)))
          end do
          do r=1,npes
            bf%a_pes(r) = bf%a_pes(r) / sumamps
          end do
        end if
      else
        bf%a_pes(in_pes) = (1.0d0,0.0d0)
      end if
      do r=1,npes
        bf%d_pes(r) = bf%a_pes(r)
      end do
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
        call propstep(swrmbf,dt,dtnext,dtdone,t,genflg,timestrt,x-stepback,reps)
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
        call propstep(swrmbf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps)
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
          call propstep(tmpbf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps)
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
          call propstep(tmpbf,dt,dtnext,dtdone,t,genflg,timestrt,x,reps)
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

  subroutine genD_big(bs, mup, muq, restart)   !   Level 1 Subroutine

    implicit none
    type(basisfn),dimension(:),intent(inout)::bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    complex(kind=8),dimension(:,:),allocatable::ovrlp_mat
    complex(kind=8),dimension(:), allocatable::zinit, zinit2, zpq
    complex(kind=8),dimension(:), allocatable:: C_k, D
    real(kind=8) :: x,y,sumamps
    integer, intent(inout) :: restart
    integer::k, j, r, ierr

    if (errorflag .ne. 0) return
    ierr = 0

    allocate(D(size(bs)), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of D_k array in genD_big"
      errorflag=1
      return
    end if 
    
    if ((basis.ne."SWTRN")) then
      do j=1,size(bs)
        if (qss==1) then
          if (npes==2) then
!            if (randfunc.eq."GAUS") then
              call random_number(x)
              call random_number(y)
              x = x * 2.0d0 * pirl
              y = y * 2.0d0 * pirl
!            else
!              x = ZBQLUAB(0.0d0,2.0d0*pirl)
!              y = ZBQLUAB(0.0d0,2.0d0*pirl)
!            end if
            do r=1,npes
              if (r==in_pes) then
                bs(j)%a_pes(r)=cmplx(dcos(x),0.0d0,kind=8)
              else
                bs(j)%a_pes(r)=cmplx(dsin(x)*dcos(y),dsin(x)*dsin(y),kind=8)
              end if
            end do
          else
            sumamps = 0.0d0
            do r=1,npes
!              if (randfunc.eq."GAUS") then
                call random_number(x)
                call random_number(y)
                x = x * 2.0d0 - 1.0d0
                y = y * 2.0d0 - 1.0d0
!              else
!                x = ZBQLUAB(-1.0d0,1.0d0)
!                y = ZBQLUAB(-1.0d0,1.0d0)
!              end if
              bs(j)%a_pes(r) = cmplx(x,y,kind=8)
              sumamps = sumamps + dble(bs(j)%a_pes(r) * dconjg(bs(j)%a_pes(r)))
            end do
            do r=1,npes
              bs(j)%a_pes(r) = bs(j)%a_pes(r) / dsqrt(sumamps)
            end do
          end if
        else
          bs(j)%a_pes(in_pes) = (1.0d0,0.0d0)
        end if
        do r=1,npes
          bs(j)%d_pes(r) = bs(j)%a_pes(r)
        end do
      end do
    end if

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
      C_k(k) = ovrlpij(zpq, zinit)
    end do
       
    do k=1,size(bs)
      C_k(k) = C_k(k) * dconjg(bs(k)%a_pes(in_pes))
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

    do k=1,size(bs)
      bs(k)%D_big = D(k)
    end do

    if (method == "MCEv1") then
      do j=1,size(bs)
        do r=1,npes
          bs(j)%a_pes(r) = bs(j)%a_pes(r) * bs(j)%D_big
          bs(j)%d_pes(r) = bs(j)%a_pes(r)
        end do
        bs(j)%D_big = (1.0d0,0.0d0)
      end do
    end if

    return

  end subroutine genD_big

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
