MODULE derivsMCE

  use globvars
  use Ham
  use alarrays
  use outputs
  use redirect
  use readpars

!***********************************************************************************!
!*
!*         MCE Derivatives Module
!*           
!*   Contains subroutines for:
!*
!*      1) Calling the time derivative subroutines for each variable and returning
!*        the results to the timestep subroutine
!*      2) Calculating the time derivative of the z values
!*      3) Calculating the derivative of the single configuration MCEv1 d prefactor
!*      4) Calculating the derivative of the single configuration MCEv2 d prefactor
!*      5) Calculating the time derivative of the classical action
!*      6) Calculating the derivative of the multi-configuration MCEv1 D prefactor
!*      7) Calculating the derivative of the multi-configuration MCEv2D prefactor
!*      8) Calculating the first component of the Hamiltonian used in 7)
!*      9) Calculating the second component of the Hamiltonian used in 7)
!*     10) Calculating the third component of the Hamiltonian used in 7) and in 3)
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

  subroutine deriv(bsin, bsout, rkstp, time, genflg, reps, x, map_bfs)

    implicit none
    type(basisfn), dimension (:), intent (in) :: bsin
    type(basisfn), dimension (:), intent (inout) :: bsout
    real(kind=8), intent(inout) :: time
    integer, intent(inout), dimension(:,:) :: map_bfs
    integer, intent(in) :: rkstp, genflg, reps, x
    
    type(basisfn), dimension(:), allocatable :: dummybs
    complex(kind=8), dimension(:,:), allocatable :: dz, dd
    complex(kind=8), dimension(:), allocatable ::dD_big
    real(kind=8), dimension(:,:), allocatable :: ds
    real(kind=8), dimension(:), allocatable :: dummy_arr
    integer :: k, m, r, ierr

    if (errorflag .ne. 0) return

    ierr = 0
    
    allocate(dummy_arr(ndim), stat = ierr)
    if (ierr/=0) then
      write (0,"(a,i0)") "Error allocating the dummy array for AIMC-MCE1 outputs. ierr was ", ierr
      errorflag = 1
      return
    end if
    dummy_arr = 0.0d0

    allocate (dz(size(bsin),ndim), stat=ierr)
    if (ierr==0) allocate (dd(size(bsin),npes),stat=ierr)
    if (ierr==0) allocate (ds(size(bsin),npes),stat=ierr)
    if (ierr==0) allocate (dD_big(size(bsin)),stat=ierr)
    if (ierr/=0) then
      write(0,"(a,i0)") "Error in derivatives array allocation in deriv. ierr was ", ierr
      errorflag=1
      return
    end if

    select case (method)
    
      case ("MCEv1")
        dz=zdot(bsin,time)
        ds=sdot(bsin,dz,time)
        dd=ddotv1(bsin,time,dz)
        dD_big=(0.0d0,0.0d0)
        if (errorflag .ne. 0) return     
        do k = 1,size(bsin)
          do m = 1,ndim
            bsout(k)%z(m) = dz(k,m)
          end do
          do r = 1,npes
            bsout(k)%d_pes(r) = dd(k,r)
            bsout(k)%s_pes(r) = ds(k,r)
          end do
          bsout(k)%D_big = dD_big(k)
        end do 
              
      case ("MCEv2")
        dz=zdot(bsin,time)
        ds=sdot(bsin,dz,time)
        dd=ddotv2(bsin,time)
        if (((basis=="TRAIN").or.(basis=="SWTRN")).and.(genflg==1)) then
          dD_big=(0.0d0,0.0d0)
        else
          dD_big=bigDdotv2 (bsin,rkstp,time,x,reps,dz)
        end if
        if (errorflag .ne. 0) return     
        do k = 1,size(bsin)
          do m = 1,ndim
            bsout(k)%z(m) = dz(k,m)
          end do
          do r = 1,npes
            bsout(k)%d_pes(r) = dd(k,r)
            bsout(k)%s_pes(r) = ds(k,r)
          end do
          bsout(k)%D_big = dD_big(k)
        end do  
        
      case ("CCS")
        dz=zdot(bsin,time)
        ds=sdot(bsin,dz,time)
        dd=(0.0d0,0.0d0)
        if (((basis=="TRAIN").or.(basis=="SWTRN")).and.(genflg==1)) then
          dD_big=(0.0d0,0.0d0)
        else
          dD_big=bigDdotv2 (bsin,rkstp,time,x,reps,dz)
        end if
        if (errorflag .ne. 0) return     
        do k = 1,size(bsin)
          do m = 1,ndim
            bsout(k)%z(m) = dz(k,m)
          end do
          do r = 1,npes
            bsout(k)%d_pes(r) = dd(k,r)
            bsout(k)%s_pes(r) = ds(k,r)
          end do
          bsout(k)%D_big = dD_big(k)
        end do
        
      case ("AIMC1")
        dz=zdot(bsin,time)
        ds=sdot(bsin,dz,time)
        dd=ddotv2(bsin,time)
        dD_big=(0.0d0,0.0d0)
        if (errorflag .ne. 0) return     
        do k = 1,size(bsin)
          do m = 1,ndim
            bsout(k)%z(m) = dz(k,m)
          end do
          do r = 1,npes
            bsout(k)%d_pes(r) = dd(k,r)
            bsout(k)%s_pes(r) = ds(k,r)
          end do
          bsout(k)%D_big = dD_big(k)
        end do
        call outbs (bsout, reps, dummy_arr, dummy_arr, time, x, rkstp)
        
      case ("AIMC2")
        call constrtrain (dummybs, x, time, reps, dummy_arr, dummy_arr ,rkstp, genflg, size(bsin), map_bfs)
        if (size(dummybs).ne.(size(bsin))) then
          write (0,"(a)") "Size of basis set in memory does not match basis set in input file"
          write (0,"(a)") "This could be caused by an incorrect handling of a cloning event"
          write (0,"(2(a,i0))") "Expected ", size(bsin), " but got ", size(dummybs)
          errorflag = 1
          return
        end if
        dD_big=bigDdotv2 (bsin,rkstp,time,x,reps,dz)
        if (errorflag .ne. 0) return     
        do k = 1,size(bsin)
          do m = 1,ndim
            bsout(k)%z(m) = dummybs(k)%z(m)
          end do
          do r = 1,npes
            bsout(k)%d_pes(r) = dummybs(k)%d_pes(r)
            bsout(k)%s_pes(r) = dummybs(k)%s_pes(r)
          end do
          bsout(k)%D_big = dD_big(k)
        end do  
        call deallocbs(dummybs)  
            
      case default  
        write(0,"(a)") "Error! Method unrecognised!"
        write(0,"(a)") "How did you even get this far?"
        errorflag = 1
        return
        
    end select                 

    deallocate (dz,dd,ds,dD_big, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in derivatives array deallocation in deriv"
      errorflag=1
      return
    end if

    return

  end subroutine deriv    

!------------------------------------------------------------------------------------

  function zdot(bsin,t)

    implicit none

    complex(kind=8), dimension (:) ,allocatable:: zdottemp
    complex(kind=8), dimension (:,:,:,:),allocatable :: dhdz
    complex(kind=8), dimension (:,:) ,allocatable:: z
    type(basisfn), dimension (:), intent (in) :: bsin
    complex(kind=8), dimension(size(bsin),ndim) :: zdot
    complex(kind=8), dimension(:),allocatable :: a, ac
    complex(kind=8) :: asum
    real(kind=8), intent (in) :: t
    integer :: k, m, r, s, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(zdottemp(ndim), stat = ierr)
    if (ierr==0) allocate (z(size(bsin),ndim), stat = ierr)
    if (ierr==0) allocate (dhdz(size(bsin),npes,npes,ndim), stat = ierr)
    if (ierr==0) allocate (a(npes), stat = ierr)
    if (ierr==0) allocate (ac(npes), stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in array allocation in zdot"
      errorflag=1
      return
    end if
    
    do k=1,size(bsin)
      do m=1,ndim
        z(k,m) = bsin(k)%z(m)
      end do
    end do
    
    call dh_dz(dhdz, z, t)

    do k=1,size(bsin)
      zdottemp = (0.0d0, 0.0d0)
      do r=1,npes
        a(r) = bsin(k)%a_pes(r)
        ac(r) = dconjg(a(r))
      end do
      asum = (0.0d0, 0.0d0)
      do r=1,npes
        asum = asum + (ac(r)*a(r))
        do s=1,npes
          do m=1,ndim
            zdottemp(m) = zdottemp(m) + (dhdz(k,r,s,m)*ac(r)*a(s))
          end do
        end do
      end do
      do m=1,ndim
        zdot(k,m)=(-1.0d0*i*zdottemp(m))/asum
      end do
    end do

    deallocate(zdottemp,z,dhdz,a,ac, stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in array deallocation in zdot"
      errorflag=1
      return
    end if

    return

  end function zdot

!------------------------------------------------------------------------------------

  function ddotv1(bsin,time,dz)

    implicit none

    type(basisfn), dimension (:), intent (in) :: bsin
    type(hamiltonian), dimension (:,:), allocatable :: H
    complex(kind=8), dimension (:,:),allocatable :: ovrlp, ovrlpin, d2H_temp
    complex(kind=8), dimension (:,:),allocatable :: zczd, ddot_temppes
    complex(kind=8), dimension (:,:,:),allocatable :: d2H
    complex(kind=8), dimension (:),allocatable :: atemp, ddot_temp, ddot_out, ddot_in
    complex(kind=8), dimension(size(bsin),npes) :: ddotv1
    complex(kind=8), dimension(:,:), intent(in) :: dz
    real(kind=8), intent (in) :: time
    integer :: j,k,r,s,nbf,ierr

    if (errorflag .ne. 0) return

    ierr = 0

    nbf = size(bsin)

    allocate(ovrlp(nbf,nbf),stat=ierr)
    if (ierr==0) allocate (ovrlpin(nbf,nbf),stat=ierr)
    if (ierr==0) allocate (zczd(nbf,nbf),stat=ierr)
    if (ierr==0) allocate (d2H_temp(nbf,nbf),stat=ierr)
    if (ierr==0) allocate (d2H(nbf,nbf,npes),stat=ierr)
    if (ierr==0) allocate (atemp(nbf),stat=ierr)
    if (ierr==0) allocate (ddot_temp(nbf),stat=ierr)
    if (ierr==0) allocate (ddot_out(nbf),stat=ierr)
    if (ierr==0) allocate (ddot_in(nbf),stat=ierr)
    if (ierr==0) allocate (ddot_temppes(nbf, npes),stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of matrices or arrays in ddotv1"
      errorflag=1
      return
    end if

    ovrlp = ovrlpmat(bsin)

    call allocham(H,nbf)
    call Hord(bsin,H,time)

    zczd = zczdot(bsin, dz)

    do r=1,npes
      do j=1,nbf
        do k=1,nbf
          if (j==k) then
            d2H(k,j,r) = (0.0d0,0.0d0)
          else
            d2H(k,j,r) = H(k,j)%Hjk(r,r) - H(j,j)%Hjk(r,r) - zczd(k,j)
          end if
        end do
      end do
    end do

    do r=1,npes
      do j=1,nbf
        do k=1,nbf
          if (k.ne.j) then
            d2H_temp(k,j) = d2H(k,j,r) * ovrlp(k,j)
          else
            d2H_temp(k,j) = (0.0d0,0.0d0)
          end if
        end do
        atemp(j) = bsin(j)%A_pes(r)
      end do
      ddot_temp = matmul(d2H_temp,atemp)
      do k=1,nbf
        ddot_temppes(k,r) = ddot_temp(k)
      end do
    end do
   

    do r=1,npes
      ddot_temp = (0.0d0,0.0d0)
      do s=1,npes
        if (s.ne.r) then
          do k=1,nbf
            do j=1,nbf
              ddot_temp(k) = ddot_temp(k)+H(k,j)%Hjk(r,s)*ovrlp(k,j)*&
                              bsin(j)%a_pes(s)
            end do
          end do
        end if
      end do
      do k=1,nbf
        ddot_temppes(k,r) = ddot_temppes(k,r) + ddot_temp(k)
      end do
    end do

    call deallocham(H)
 
    do r=1,npes   

      do k=1,nbf
        ddot_in(k) = ddot_temppes(k,r)
        do j=1,nbf
          ovrlpin(j,k) = ovrlp(j,k)    !!!! ovrlp is changed by the matinv subroutine
        end do                         !!!! so must be redefined for each pes
      end do

      if (matfun.eq.'zgesv') then
        call lineq (ovrlpin, ddot_in, ddot_out)
      else if (matfun.eq.'zheev') then
        call matinv2(ovrlpin, ddot_in, ddot_out)
      else
        write(0,"(a)") "Error! Matrix function not recognised! Value is ", matfun
        errorflag = 1
        return
      end if 

      do k=1,nbf
        ddotv1(k,r)=-1.0d0*i*ddot_out(k)*cdexp(-1.0d0*i*(bsin(k)%s_pes(r)))
      end do
    
    end do

    deallocate(ovrlp,ovrlpin,zczd,d2H_temp,d2H,atemp,stat=ierr)
    if(ierr==0) deallocate(ddot_temp,ddot_out,ddot_in,ddot_temppes,stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of matrices or arrays in ddotv1"
      errorflag=1
      return
    end if
 
    return
 
  end function ddotv1  

!------------------------------------------------------------------------------------

  function ddotv2(bsin,t)

    implicit none

    type(basisfn), dimension (:), intent (in) :: bsin
    real(kind=8), intent (in) :: t
    complex(kind=8), dimension(:,:,:),allocatable :: Hkk
    complex(kind=8), dimension(:),allocatable :: dk
    complex(kind=8), dimension(size(bsin),npes) :: ddotv2
    complex(kind=8), dimension(:,:),allocatable :: z
    real(kind=8), dimension(:),allocatable :: Sk
    integer :: k, m, r, s, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    allocate(Hkk(size(bsin),npes,npes),dk(npes),z(size(bsin),ndim),Sk(npes),stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of matrices or arrays in ddotv2"
      errorflag=1
      return
    end if

    do r=1,size(ddotv2,2)
      do k=1,size(ddotv2,1)
        ddotv2(k,r) = (0.0d0, 0.0d0)
      end do
    end do

    do k=1,size(bsin)
      do m=1,ndim
        z(k,m) = bsin(k)%z(m)
      end do
    end do
      
    call Hijdiag(Hkk,z,t)
    
    do k=1,size(bsin)
      dk = bsin(k)%d_pes
      Sk = bsin(k)%S_pes
      do r=1,npes
        do s=1,npes
          if (r.ne.s) then
            ddotv2(k,r) = ddotv2(k,r) + Hkk(k,r,s) * dk(s) &!* ovrlpij(z,z) &
                                                 * cdexp(i*(Sk(s)-Sk(r)))
          end if
        end do
      end do
    end do

    do r=1,size(ddotv2,2)
      do k=1,size(ddotv2,1)
        ddotv2(k,r) = -1.0d0*i*ddotv2(k,r)
      end do
    end do

    deallocate(Hkk,dk,z,Sk,stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of matrices or arrays in ddotv2"
      errorflag=1
      return
    end if

    return

  end function ddotv2

!------------------------------------------------------------------------------------

  function sdot(bsin,dz,t)

    implicit none

    type(basisfn), dimension (:), intent (in) :: bsin
    complex(kind=8), dimension (:,:), intent(in) :: dz
    real(kind=8), intent (in) :: t
    real(kind=8), dimension (size(bsin),npes) :: sdot
    complex(kind=8), dimension(:,:), allocatable :: zk, zkc, zkdot, zkdotc
    complex(kind=8), dimension(:,:,:), allocatable :: Hkk
    integer :: k, r, m, ierr
    complex(kind=8) :: zsum

    if (errorflag .ne. 0) return
    
    ierr = 0

    allocate(zk(size(bsin),ndim), stat=ierr)
    if (ierr==0) allocate(zkc(size(bsin),ndim), stat=ierr)
    if (ierr==0) allocate(zkdot(size(bsin),ndim), stat=ierr)
    if (ierr==0) allocate(zkdotc(size(bsin),ndim), stat=ierr)
    if (ierr==0) allocate(Hkk(size(bsin),npes,npes),stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of matrices or arrays in sdot"
      errorflag=1
      return
    end if

    do k = 1,size(bsin)
      do m = 1,ndim
        zk(k,m) = bsin(k)%z(m)
        zkdot(k,m) = dz(k,m)
        zkc(k,m) = dconjg(zk(k,m))
        zkdotc(k,m) = dconjg(zkdot(k,m))
      end do
    end do
    
    call Hijdiag(Hkk,zk,t)
      
    do k=1,size(bsin)
      zsum = (0.0d0, 0.0d0)
      do m = 1,ndim
        zsum = zsum + i*0.5d0*(zkdot(k,m)*zkc(k,m)-zkdotc(k,m)*zk(k,m))
      end do
      do r = 1,npes
        sdot(k,r) = zsum - dble(Hkk(k,r,r))
      end do
    end do

    deallocate(zk,zkc,zkdot,zkdotc,Hkk,stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of matrices or arrays in sdot"
      errorflag=1
      return
    end if

    return

  end function sdot

!------------------------------------------------------------------------------------

  function bigDdotv1(nbf)
 
    implicit none

    integer, intent (in) :: nbf
    complex(kind=8), dimension(nbf)::bigDdotv1
    integer :: j

    if (errorflag .ne. 0) return

    do j=1,nbf
      bigDdotv1 = (0.0d0,0.0d0)
    end do

    return

  end function bigDdotv1    

!------------------------------------------------------------------------------------

  function bigDdotv2(bsin,rkstp, time,x,reps,dz)

    implicit none

    type(basisfn), dimension (:), intent (in) :: bsin 
    type(hamiltonian), dimension (:,:), allocatable :: H
    complex(kind=8), dimension (:,:), allocatable :: ovrlp, ovrlpphi
    complex(kind=8), dimension (:,:), allocatable :: h_av_jj, h_av_jk, zconjzdot
    complex(kind=8), dimension (:,:), allocatable :: d2H, d2H2, d2Hdiff
    complex(kind=8), dimension(:), allocatable::tempD, chk, Ddot, d2HD, Dtemp
    complex(kind=8), dimension(size(bsin)) :: bigDdotv2
    complex(kind=8), dimension(:,:), intent(in) :: dz
    complex(kind=8) :: ovrlpdif, asum
    real(kind=8), intent(in) :: time
    integer, intent(in) :: rkstp, x, reps
    integer :: j, k, nbf, ierr, r

    if (errorflag .ne. 0) return

    ierr = 0

    nbf = size(bsin)

    allocate (ovrlp(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (ovrlpphi(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (h_av_jj(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (h_av_jk(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (zconjzdot(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (d2H(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (d2H2(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (d2Hdiff(nbf,nbf),stat=ierr)
    if (ierr == 0) allocate (tempD(nbf),stat=ierr)
    if (ierr == 0) allocate (chk(nbf),stat=ierr)
    if (ierr == 0) allocate (Ddot(nbf),stat=ierr)
    if (ierr == 0) allocate (d2HD(nbf),stat=ierr)
    if (ierr == 0) allocate (Dtemp(nbf),stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of matrices or arrays in bigDdotv2"
      errorflag=1
      return
    end if    

    ovrlp = ovrlpmat(bsin)
    do k=1,nbf
      do j=k,nbf
        asum = (0.0d0,0.0d0)
        do r=1,npes
          asum = asum + (dconjg(bsin(j)%a_pes(r))*bsin(k)%a_pes(r))
        end do
        ovrlpphi(j,k)=ovrlp(j,k)*asum
        if (j.ne.k) then
          ovrlpphi(k,j)=dconjg(ovrlpphi(j,k))
        end if
      end do
    end do
    
    if ((rkstp.eq.1).and.(time.eq.0.0d0)) then
      do j=1,nbf
        do k=1,nbf
          ovrlpdif = ovrlpphi(j,k) - ovrlp(j,k) 
          if ((ovrlpdif.ne.0.0d0).and.(basis.ne."TRAIN").and.(basis.ne."SWTRN").and.&
            & (cloneflg.ne."BLIND").and.(cloneflg.ne."BLIND+")) then
            write(0,"(a)") "Error! Initial phi-overlap has disimilarilies to z-overlap"
            write(0,'(a,a,i0,a,i0,a)'), "These matrices should be identical ",&
                                  "but differences found at coordinate ", j,",",k,"."
            write(0,'(a,4(e15.8,a))'), "Expected (",dimag(i*ovrlp(j,k)),","&
                                 ,dimag(ovrlp(j,k)),") but got (",&
                                  dimag(i*ovrlpphi(j,k)),",",dimag(ovrlpphi(j,k)),")"
            errorflag = 1
          end if
        end do
      end do
      if (errorflag == 1) return
    end if

    call allocham(H, nbf)
    call Hord(bsin, H, time)
    
    h_av_jj = Hjk_avrg(H,bsin)
    h_av_jk = phiHphi(H,bsin)
    zconjzdot = zczdot(bsin,dz)
    if (errorflag .ne. 0) return    

    do k=1,nbf
      Dtemp(k) = bsin(k)%D_big
      do j=1,nbf
        d2H(j,k) = h_av_jk(j,k) - h_av_jj(j,k) - zconjzdot(j,k)
        if (d2H(j,k)/=d2H(j,k)) then
          write(0,"(a,i4,a,i4,a)") "d2H(", j, ",", k, ") is NaN. Terminating"
          errorflag = 1
          return
        end if
      end do
    end do
    
    if ((rkstp.eq.1).and.(mod((x-1),100)==0).and.(basis.ne."TRAIN").and.(sys.eq."SB").and.(debug==1)) then
      call outd2hovrlp(d2H,ovrlp,x,reps)
    end if
    
    do k=1,nbf
      do j=1,nbf
        d2H(j,k) = d2H(j,k) * ovrlp(j,k)
      end do
    end do

    tempD = matmul(d2H, Dtemp)

    do k=1,nbf
      tempD(k)=-1.0d0*i*tempD(k)
    end do


    do j=1,nbf
      Ddot(j) = (0.0d0, 0.0d0)
      if (tempD(j)/=tempD(j)) then
        write(0,"(a,i4,a)") "tempD(", j, ") is NaN. Terminating"
        errorflag = 1
      return
      end if
    end do

    d2hD = tempD

    if (matfun.eq.'zgesv') then
      call lineq (ovrlpphi, tempD, Ddot)
    else if (matfun.eq.'zheev') then
      call matinv2(ovrlpphi, tempD, Ddot)
    else
      write(0,"(a)") "Error! Matrix function not recognised! Value is ", matfun
      errorflag = 1
      return
    end if

    do k=1,nbf
      bigDdotv2(k)=Ddot(k)
    end do

   call deallocham(H)

    deallocate (ovrlp,ovrlpphi,h_av_jj,h_av_jk,zconjzdot,d2H,d2H2,d2Hdiff,stat=ierr)
    if (ierr == 0 ) deallocate (tempD,chk,Ddot,d2HD,Dtemp,stat = ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in deallocation of matrices or arrays in bigDdotv2"
      errorflag=1
      return
    end if   

    return

  end function bigDdotv2

!------------------------------------------------------------------------------------

  function phiHphi(H,bsin)

    implicit none
    type(basisfn), dimension (:), intent (in) :: bsin
    type(hamiltonian), dimension (:,:), intent(in) :: H
    complex(kind=8), dimension (size(bsin),size(bsin)) :: phiHphi
    integer :: j, k , r, s

    if (errorflag .ne. 0) return

    do k=1,size(phiHphi,2)
      do j=1,size(phiHphi,1)
        phiHphi(j,k) = (0.0d0, 0.0d0)
        do r=1,npes
          do s=1,npes
            phiHphi(j,k) = phiHphi(j,k) + (dconjg(bsin(j)%a_pes(r)) * &
                                        H(j,k)%Hjk(r,s) * bsin(k)%a_pes(s))
          end do
        end do
        if (phiHphi(j,k)/=phiHphi(j,k)) then
          write(0,"(a,i4,a,i4,a)") "phiHphi(", j, ",", k, ") is NaN. Terminating"
          errorflag = 1
          return
        end if
      end do
    end do

    return

  end function phiHphi

!------------------------------------------------------------------------------------

  function Hjk_avrg(H,bsin)

    implicit none
    type(basisfn), dimension (:), intent (in) :: bsin
    type(hamiltonian), dimension (:,:), intent(in) :: H
    complex(kind=8), dimension (size(bsin),size(bsin)) :: Hjk_avrg
    integer :: j, k , r, s

    if (errorflag .ne. 0) return

    do k=1,size(Hjk_avrg,2)
      do j=1,size(Hjk_avrg,1)
        Hjk_avrg(j,k) = (0.0d0, 0.0d0)
        do r=1,npes
          do s=1,npes
            Hjk_avrg(j,k) = Hjk_avrg(j,k) + (dconjg(bsin(j)%a_pes(r)) * bsin(k)%a_pes(s) *&
                                              H(k,k)%Hjk(r,s))
          end do
        end do
        if (Hjk_avrg(j,k)/=Hjk_avrg(j,k)) then
          write(0,"(a,i4,a,i4,a)") "Hjk_avrg(", j, ",", k, ") is NaN. Terminating"
          errorflag = 1
          return
        end if
      end do
    end do

    return

  end function Hjk_avrg

!------------------------------------------------------------------------------------

  function zczdot(bsin,dz)

    type(basisfn), dimension (:), intent (in) :: bsin 
    complex(kind=8), dimension (size(bsin),size(bsin)) :: zczdot
    complex(kind=8), dimension (:,:), intent(in) :: dz
    complex(kind=8) :: asum
    integer :: j, k, m, r, ierr

    if (errorflag .ne. 0) return

    ierr = 0

    do k=1,size(zczdot,2)
      do j=1,size(zczdot,1)
        zczdot(j,k) = (0.0d0, 0.0d0)
        do m=1,ndim
          zczdot(j,k) = zczdot(j,k) + ((dconjg(bsin(j)%z(m)-bsin(k)%z(m)))*dz(k,m))
        end do
        asum = (0.0d0,0.0d0)
        do r=1,npes
          asum = asum + (dconjg(bsin(j)%a_pes(r)) * bsin(k)%a_pes(r))
        end do                                        
        if (zczdot(j,k)/=zczdot(j,k)) then
          write(0,"(a,i4,a,i4,a)") "zczdot(", j, ",", k, ") is NaN. Terminating"
          errorflag = 1
          return
        end if
        if (method == "MCEv1") then
          zczdot(j,k) = i * zczdot(j,k)
        else
          zczdot(j,k) = i * zczdot(j,k) * asum
        end if
      end do
    end do 

    return

  end function zczdot

!***********************************************************************************!

end module derivsMCE
