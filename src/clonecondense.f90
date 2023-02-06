module clonecondense
  use alarrays
  use outputs
  use readpars
  use Ham

!*****************************************************************
!*   First attempt at translating the python file clonecombine.py into fortran to be called automatically
!*   and interact with bases much easier, as well as using pre-written functions.
!*
!*    Contains subroutines for:
!*    
!*     1) Recondensing the clones that are created through version 1 propagation
!*    
!*
!*****************************************************************

contains 
  

  subroutine clone_condense(bsetarr,dt,tsstart,tsend,reps,nclones,nbf,absnorm,acf_t,extra,absehr,pops,mup,muq)
    implicit none

    type(basisset), dimension(:,:), intent(inout) :: bsetarr
    real(kind=8), dimension(:), intent(inout)  :: absnorm, absehr
    complex(kind=8), dimension (:), intent(inout):: acf_t, extra
    real(kind=8), dimension (:,:), intent(inout) :: pops
    real(kind=8), dimension(:,:), allocatable :: populations, ctarray, normpfs
    real(kind=8) :: pophold1, pophold2, ehrtmp
    integer(kind=4), intent(in):: tsstart, tsend, nclones, reps, nbf
    real(kind=8), dimension(:), intent(inout) :: mup,muq
    real(kind=8)::crossterm1, crossterm2, time, rescale, nrmtmp
    real(kind=8), intent(in):: dt
    integer:: h, f, i, e, t, range, j, n
    complex(kind=8), dimension(:,:), allocatable :: ovrlp
    character(LEN=4):: rep
    complex(kind=8)::normtemp,acft, ehren, extmp 
    
    ! need to have a do loop for each time a bset is passed between the array and an actionable basis set.?????

    ! Defining value of varibles
    crossterm1 = 0
    crossterm2 = 0
    ! time = dt * tsstart
    range = tsend-tsstart
  
    allocate(normpfs(13,range))
    allocate(populations(range,2))
    allocate(ctarray(range,2))
    allocate(ovrlp(nbf,nbf))  
    ovrlp = 0.0
    populations = 0.0
    ctarray = 0.0
    normpfs = 0.0
    
    time = dt*tsstart+dt
    do t=1,range
      ehren = (0.0d0,0.0d0)
      do i=1, nclones
        ovrlp=ovrlpmat(bsetarr(i,t)%bs)
        pophold1 = pop(bsetarr(i,t)%bs,1,ovrlp)
        pophold2 = pop(bsetarr(i,t)%bs,2,ovrlp)
        normtemp = norm(bsetarr(i,t)%bs,ovrlp)
        acft = acf(bsetarr(i,t)%bs,mup,muq)
        acf_t(t+tsstart) = acf_t(t+tsstart) + acft
        call extras(extmp, bsetarr(i,t)%bs)
        extra(t+tsstart) = extra(t+tsstart) + extmp
        do j = 1,nbf
          ehren = ehren + HEhr(bsetarr(i,t)%bs(j), time, reps)
        end do
        ehrtmp = abs(ehren)
        absehr(t+tsstart) = absehr(t+tsstart) + ehrtmp
        nrmtmp = 1/nclones !sqrt(dble(normtemp*dconjg(normtemp)))
        absnorm(t+tsstart) = absnorm(t+tsstart) + nrmtmp
        populations(t,1) = populations(t,1) + pophold1
        populations(t,2) = populations(t,2) + pophold2
        !write(6,*) 'at time', t, 'populations are', pop(bsetarr(i,t)%bs,1,ovrlp),pop(bsetarr(i,t)%bs,2,ovrlp)
        do j=i+1,nclones
          call cross_terms(bsetarr(i,t)%bs,bsetarr(j,t)%bs, crossterm1,crossterm2,nbf)
          ctarray(t,1) = ctarray(t,1) + crossterm1
          ctarray(t,2) = ctarray(t,2) + crossterm2
        end do 
      end do 
      
      normpfs(1,t) = time
      normpfs(2,t) = 1
      normpfs(3,t) = dble(acf_t(t+tsstart))
      normpfs(4,t) = dimag(acf_t(t+tsstart))
      normpfs(5,t) = abs(acf_t(t+tsstart))
      normpfs(6,t) = dble(extra(t+tsstart)) 
      normpfs(7,t) = dimag(extra(t+tsstart))
      normpfs(8,t) = abs(extra(t+tsstart))
      normpfs(9,t) = absehr(t+tsstart)
      normpfs(10,t) = populations(t,1) + ctarray(t,1)
      normpfs(11,t) = populations(t,2) + ctarray(t,2)
      rescale = normpfs(10,t) + normpfs(11,t)
      normpfs(10,t) = normpfs(10,t)/rescale
      normpfs(11,t) = normpfs(11,t)/rescale
      open(21061,file='popfromcc.out',access='append')
      write(21061,*) populations(t,1),populations(t,2), ctarray(t,1), ctarray(t,2), rescale, normpfs(10,t), normpfs(11,t)
      close(21061)
      pops(t+tsstart+1,1) = pops(t+tsstart+1,1) + normpfs(10,t)
      pops(t+tsstart+1,2) = normpfs(11,t) + pops(t+tsstart+1,2)
      normpfs(12,t) = normpfs(10,t) + normpfs(11,t)
      normpfs(13,t) = normpfs(10,t) - normpfs(11,t)
      time = time + dt
    end do    
    normpfs(2,1) = 1
    absnorm(tsstart+1) = absnorm(tsstart+1)/2
    ! Now to put that the population matrix into a collated output file.
    write(rep,"(i4.4)") reps
    open(11204, file = "normpop-"//trim(rep)//".out",access='append')
    do n = 1, range
     write (11204,"(13(1x,es16.8e3))") normpfs(:,n)
    end do 
    close(11204)


    write(6,*) 'CONDENSING FINISHED'
    
  end subroutine

  subroutine cross_terms(bs1,bs2, crossterm1, crossterm2,nbf)
    type(basisfn), dimension (:), intent(in):: bs1, bs2
    real(kind=8), intent(inout) ::  crossterm1, crossterm2
    integer, intent(in) :: nbf
    integer :: j, k, ierr
    complex(kind=8), dimension(nbf,nbf) ::ovrlp12, ovrlp21
    complex(kind=8) :: sumamps1, sumamps2, ovrlp1, ovrlp2, cpop11, cpop12, cpop21,cpop22
    



    sumamps1 = (0.0d0, 0.0d0)
    sumamps2 = (0.0d0, 0.0d0)
    crossterm1 = 0
    crossterm2 = 0
    ovrlp1 = (0.0d0,0.0d0)
    ovrlp2 = (0.d0,0.0d0)
    cpop11 = (0.d0,0.0d0)
    cpop12 = (0.d0,0.0d0)
    cpop21 = (0.d0,0.0d0)
    cpop22 = (0.d0,0.0d0)
  
    !Calculating the overlap between the basis sets and summing over amplitudes
    
    ovrlp12 = ovrlpmatcross(bs1,bs2)
    ovrlp21=transpose(conjg(ovrlp12)) 

    ! Calculating the cross terms

    do k=1,nbf
      do j=1,nbf
        cpop11 = cpop11 + dconjg(bs1(j)%D_big) * ovrlp12(j,k) * bs2(k)%D_big &
                * dconjg(bs1(j)%d_pes(1)) * bs2(k)%d_pes(1) &
                * cdexp(i*(bs1(k)%s_pes(1)-bs2(j)%s_pes(1)))
      end do
    end do
    do k=1,nbf
      do j=1,nbf
        cpop12 = cpop12 + dconjg(bs2(j)%D_big) * ovrlp21(j,k) * bs1(k)%D_big &
                * dconjg(bs2(j)%d_pes(1)) * bs1(k)%d_pes(1) &
                * cdexp(i*(bs2(k)%s_pes(1)-bs1(j)%s_pes(1)))
       
      end do
    end do
    do k=1,nbf
      do j=1,nbf
        cpop21 = cpop21 + dconjg(bs1(j)%D_big) * ovrlp12(j,k) * bs2(k)%D_big &
                * dconjg(bs1(j)%d_pes(2)) * bs2(k)%d_pes(2) &
                * cdexp(i*(bs1(k)%s_pes(2)-bs2(j)%s_pes(2)))

      end do
    end do

    do k=1,nbf
      do j=1,nbf
        cpop22 = cpop22 + dconjg(bs2(j)%D_big) * ovrlp21(j,k) * bs1(k)%D_big &
                * dconjg(bs2(j)%d_pes(2)) * bs1(k)%d_pes(2) &
                * cdexp(i*(bs2(k)%s_pes(2)-bs1(j)%s_pes(2)))

      end do
    end do

    crossterm1 = real(cpop11+cpop12) 
    crossterm2 = real(cpop21+cpop22)

    ! write(6,*) 'pops'
    ! print*,cpop11 
    ! print*,cpop12
    ! print*,cpop21 
    ! print*,cpop22 


    ! write(6,*) 'crossterms, ', crossterm1, crossterm2
   
    
  

  end subroutine cross_terms

  function ovrlpmatcross(bs1,bs2)

    implicit none
    type(basisfn), dimension (:), intent(in)::bs1, bs2
    complex(kind=8), dimension(size(bs1),size(bs1))::output,ovrlpmatcross
    integer::k,j,n

    if (errorflag .ne. 0) return

    ovrlpmatcross = (0.0d0, 0.0d0)

    n = size(bs1)

    do k=1,n
      do j=1,n 
        output(j,k)=ovrlpij(bs2(j)%z,bs1(k)%z)
      end do
    end do
    
    ovrlpmatcross=output
    
    return

  end function ovrlpmatcross

End module clonecondense