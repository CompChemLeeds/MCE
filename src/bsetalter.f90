MODULE bsetalter

  use globvars
  use Ham
  use alarrays
  use outputs

!***********************************************************************************!
!*
!*         Basis Set alteration module Module
!*
!*   Contains subroutines for:
!*
!*      1) Relocating the wavefunction on a static grid
!*      2) Leaking the basis set, used for large q values in the Hennon-Heiles model
!*      3) Cloning subroutine, which increases the basis set, used mainly for the
!*              spin boson model
!*      4) Retriving cloning information from previous calculations when using the
!*              AIMC-MCE propagation system
!*
!***********************************************************************************!

contains


!--------------------------------------------------------------------------------------------------

  subroutine reloc_basis(bsnew, bsold, x) ! Relocation for MCEv1 or MCEv2 after cloning

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bsnew   ! This is the new basis set
    type(basisfn), dimension(:), allocatable, intent(inout) :: bsold   ! This is the previous basis set.
    integer, intent(in) :: x

    complex(kind=8), dimension (:,:), allocatable :: ovrlp_mat
    complex(kind=8), dimension (:), allocatable :: cnew, dnew
    complex(kind=8) :: sumamps
    integer :: j, k, r, ierr, nbfnew, nbfold

    if (errorflag .ne. 0) return

    nbfnew = size(bsnew)
    nbfold = size(bsold)

    if (method=="MCEv1") then
      do j=1,nbfold
        sumamps = (0.0d0, 0.0d0)
        do r=1,npes
          sumamps = sumamps + dconjg(bsold(j)%a_pes(r))*bsold(j)%a_pes(r)
        end do
        bsold(j)%D_big = cdsqrt(sumamps)
        do r=1,npes
          bsold(j)%d_pes(r) = bsold(j)%d_pes(r)/bsold(j)%D_big
          bsold(j)%a_pes(r) = bsold(j)%d_pes(r) * cdexp (i*bsold(j)%s_pes(r))
        end do
      end do
      do j=1,nbfnew
        sumamps = (0.0d0, 0.0d0)
        do r=1,npes
          sumamps = sumamps + dconjg(bsnew(j)%a_pes(r))*bsnew(j)%a_pes(r)
        end do
        bsnew(j)%D_big = cdsqrt(sumamps)
        do r=1,npes
          bsnew(j)%d_pes(r) = bsnew(j)%d_pes(r)/bsnew(j)%D_big
          bsnew(j)%a_pes(r) = bsnew(j)%d_pes(r) * cdexp (i*bsnew(j)%s_pes(r))
        end do
      end do
    end if

    allocate (cnew(nbfnew), stat=ierr)
    if (ierr==0) allocate (ovrlp_mat(nbfnew,nbfold), stat=ierr)
    if (ierr/=0) then
      write (0,"(a)") "Error in allocating temporary z values for reprojection subroutine"
      errorflag=1
      return
    end if

    ! Ovrlp_mat is the overlap with the initial wavepacket
    do j=1,nbfold
      do k=1,nbfnew
        cnew(k) = (0.0d0, 0.0d0)
        sumamps = (0.0d0, 0.0d0)
        do r=1,npes
          sumamps = sumamps + dconjg(bsnew(k)%a_pes(r))*bsold(j)%a_pes(r)
        end do
        ovrlp_mat(k,j) = ovrlpij(bsnew(k)%z(:), bsold(j)%z(:)) * sumamps
      end do
    end do

    do r=1,npes
      cnew(:) = matmul(ovrlp_mat,bsold(:)%D_big)
    end do

    deallocate(ovrlp_mat, stat=ierr)
    if (ierr==0) allocate (ovrlp_mat(nbfnew,nbfnew), stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error in de- and re-allocation of first overlap matrix in reprojection subroutine"
      errorflag = 1
      return
    end if

    ovrlp_mat = ovrlpphimat(bsnew)

    allocate (dnew(nbfnew), stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error in allocation of dnew in reprojection subroutine"
      errorflag=1
      return
    end if

    do r=1,npes
      call lineq(ovrlp_mat, cnew, dnew)
    end do

    deallocate(ovrlp_mat, stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error in deallocation of input overlap"
      errorflag=1
      return
    end if

    do k=1,nbfnew
      bsnew(k)%D_big = dnew(k)
    end do

    if (method=="MCEv1") then
      do j=1,nbfnew
        do r=1,npes
          bsnew(j)%d_pes(r) = bsnew(j)%d_pes(r)*bsnew(j)%D_big
          bsnew(j)%a_pes(r) = bsnew(j)%d_pes(r) * cdexp (i*bsnew(j)%s_pes(r))
        end do
        bsnew(j)%D_big = (1.0d0,0.0d0)
      end do
      do j=1,nbfold
        do r=1,npes
          bsold(j)%d_pes(r) = bsold(j)%d_pes(r)*bsold(j)%D_big
          bsold(j)%a_pes(r) = bsold(j)%d_pes(r) * cdexp (i*bsold(j)%s_pes(r))
        end do
        bsold(j)%D_big = (1.0d0,0.0d0)
      end do
    end if

    deallocate(dnew, stat=ierr)
    if (ierr/=0) then
      write (0,"(a)"),"Error deallocating dnew arrays in reloc"
      errorflag = 1
      return
    end if

    return

  end subroutine reloc_basis

!--------------------------------------------------------------------------------------------------

  subroutine cloning(bs,nbf,x,time,clone,clonenum,reps,mup,muq)

    !!NOTE: Cloning is only set up for 2 PESs. This needs to be generalised!

    implicit none

    real(kind=8), dimension(:), intent(in) :: mup, muq

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    type(basisfn), dimension(:), allocatable :: bsnew, bsnew2
    real(kind=8), intent(in) :: time
    integer, dimension(:), allocatable, intent(inout) :: clone, clonenum
    integer, intent (inout) :: nbf
    integer, intent (in) :: x, reps
    complex(kind=8), dimension (:), allocatable :: dz
    real(kind=8), dimension (:), allocatable :: dummy_arr, length, phi
    real(kind=8) :: brforce, normar, sumamps, p, q, pqsqrd, deltaprob, phaseval, dist
    integer, dimension(:), allocatable :: clonehere, clonecopy, clonecopy2
    integer :: k, m, j, n, nbfnew, ierr, r, s, clonetype, dummy1, dummy2, v1flg,randomintV1
    character(LEN=3)::rep
    character(LEN=4)::crrntrep,ranomintV1char
  
  

    if (errorflag==1) return

    if ((cloneflg=="YES").or.(cloneflg=="QSC").or.(cloneflg=="V1").or.((cloneflg=="BLIND+").and.(x.ne.0))) then
      clonetype = 1 ! 1=conditional cloning, 2=blind cloning
      v1flg=0 !v1flg to allow v1 cloning
    else if ((cloneflg=="BLIND").or.((cloneflg=="BLIND+").and.(x.eq.0))) then
      clonetype = 2
    else
      write (0,"(2a)") "Cloneflg is invalid. Should be 'YES', 'V1', 'BLIND' or 'BLIND+', but was ", cloneflg
      errorflag = 1
      return
    end if

    allocate (clonehere(nbf), stat=ierr)
    if (ierr==0) allocate(clonecopy(nbf), stat=ierr)
    if(cloneflg.ne.("V1")) then
      if (ierr==0) allocate(clonecopy2(nbf), stat=ierr)
    end if
    if (ierr/=0) then
      write (0,"(a)") "Error allocating the clonehere array"
      errorflag = 1
      return
    end if



    if (npes.ne.2) then
      write(6,*) "Error. Cloning currently only valid for npes=2"
      errorflag = 1
      return
    end if

    if(cloneflg=="V1") then 
      if(mod(x,clonefreq)==0)then
        do k=1,nbf 
          clone(k)=0
          clonehere(k)=0
        end do
      else 
        do k=1,nbf
          clonehere(k)=0
          clone(k)=clone(k)
        end do
      end if
    else
      do k=1,nbf
        if (clone(k).lt.x-clonefreq) clone(k)=0
        clonehere(k) = 0
      end do
    end if

    ! Build the map of which trajectories to clone, done on the basis function by basis function level (clonetype=1),
    ! or the entire basis set at once (clonetype=2)

    if (clonetype==1) then
      do k=1,nbf
        normar = 0.0d0
        do r=1,npes
          normar = normar + dconjg(bs(k)%a_pes(r))*bs(k)%a_pes(r)
        end do
        !!!!! The line below needs changing to acount for multiple PESs
        brforce = ((abs(bs(k)%a_pes(1)*bs(k)%a_pes(2))**2.0)/(normar**2.0))
        if ((brforce.gt.thresh).and.(clone(k)==0).and.(clonenum(k).lt.clonemax)) then
          !print*,x
          clone(k) = x
          clonehere(k) = 1
          if(cloneflg=="V1") then
            v1flg=1 !Set to 1 meaning at least 1 basis function is to be cloned
          end if
        end if
        clonecopy(k) = clone(k)
        if(cloneflg.ne."V1") then
          clonecopy2(k) = clonenum(k)
        end if
      end do
    else if (clonetype==2) then
      if (mod(x,clonefreq)==0) then
        do k=1,nbf
          if (clonenum(k).lt.clonemax) then
            clone(k) = x
            clonehere(k) = 1
          end if
        end do
      end if
      do k=1,nbf
        clonecopy(k) = clone(k)
        clonecopy2(k) = clonenum(k)
      end do

    end if

    ! build new sized clone mapping arrays to use for next itteration
    if(cloneflg.ne."V1") then
      nbfnew = nbf + sum(clonehere(:))

      deallocate (clone, stat=ierr)
      if (ierr==0) deallocate (clonenum, stat=ierr)
      if (ierr==0) allocate (clone(nbfnew), stat=ierr)
      if (ierr==0) allocate (clonenum(nbfnew), stat=ierr)
      if (ierr/=0) then
        write (0,"(a)") "Error in de- and re-allocation of clone arrays"
        errorflag = 1
        return
      end if
      do k=1,nbf
        clone(k) = clonecopy(k)
        clonenum(k) = clonecopy2(k)
      end do
      deallocate(clonecopy, stat=ierr)
      if (ierr==0) deallocate(clonecopy2, stat=ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error deallocating the cloning copy arrays"
        errorflag = 1
        return
      end if
    else
      nbfnew = nbf
    end if 

    ! Actual cloning section

    if ((nbfnew/=nbf).or.(v1flg==1)) then
      call allocbs(bsnew, nbfnew)
      if(v1flg==1) then
        call allocbs(bsnew2,nbfnew)
      end if

      j=1

      write(rep,"(i3.3)") reps
      open(unit=47756,file="Clonetrack-"//trim(rep)//".out",status="old",access="append",iostat=ierr)

      allocate (dz(ndim), stat=ierr)
      if (ierr/=0) then
        write(0,'(a)') "Error in allocating dz array for MCEv1 cloning/relocation"
        errorflag = 1
        return
      end if

      if (cloneflg=="QSC") then
        clonetype = 3
      else if (method=="MCEv2") then
        clonetype = 1
      else if (method=="MCEv1") then
        if(cloneflg=="V1") then
          clonetype = 4
        else
          clonetype = 2
        end if
      end if

      if(clonetype.ne.4) then      
        do k=1,nbf

          if (clonehere(k) == 1) then
            clone(k) = x
            clone(nbf+j) = x
            clonenum(k) = clonenum(k) + 1
            clonenum(nbf+j) = clonenum(k)

            if (clonetype==1) then

              ! First child trajectory
              bsnew(k)%D_big = bs(k)%D_big * abs(bs(k)%a_pes(in_pes))
              bsnew(k)%d_pes(in_pes) = bs(k)%d_pes(in_pes)/abs(bs(k)%a_pes(in_pes))
              do r=1,npes
                if (r.ne.in_pes) then
                  bsnew(k)%d_pes(r) = (0.0d0,0.0d0)
                end if
                bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
                bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
              end do
              do m=1,ndim
                bsnew(k)%z(m) = bs(k)%z(m)
              end do

              ! Second child trajectory
              bsnew(nbf+j)%D_big = bs(k)%D_big * sqrt(1.-(dconjg(bs(k)%a_pes(in_pes))*bs(k)%a_pes(in_pes)))
              bsnew(nbf+j)%d_pes(in_pes) = (0.0d0,0.0d0)
              do r=1,npes
                if (r.ne.in_pes) then
                  if (x.eq.0) then
                    bsnew(nbf+j)%d_pes(r) = (1.0d0,0.0d0)
                  else
                    bsnew(nbf+j)%d_pes(r) = bs(k)%d_pes(r)/&
                                    sqrt(1.-(dconjg(bs(k)%a_pes(in_pes))*bs(k)%a_pes(in_pes)))
                  end if
                end if
                bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
                bsnew(nbf+j)%a_pes(r) = bsnew(nbf+j)%d_pes(r) * cdexp(i*bsnew(nbf+j)%s_pes(r))
              end do
              do m=1,ndim
                bsnew(nbf+j)%z(m) = bs(k)%z(m)
              end do

            else if (clonetype==2) then

              !!!!!!!!!!NOTE!!!!!!!!!!!

              !This is only half of the MCEv1 cloning procedure. After cloning occurs,
              !the reloc_basis subroutine is called to recalculate the amplitudes.
              ! Checking that the reloc function works as expected would be a good place to start

              ! First child trajectory
              bsnew(k)%D_big = bs(k)%D_big
              bsnew(k)%d_pes(in_pes) = bs(k)%d_pes(in_pes)
              do m=1,ndim
                dz(m)=cmplx(ZBQLNOR(dble(bs(k)%z(m)),sqrt(0.5d0)),ZBQLNOR(dimag(bs(k)%z(m)),sqrt(0.5)))
              end do
              do r=1,npes
                if (r.ne.in_pes) then
                  bsnew(k)%d_pes(r) = (0.0d0,0.0d0)
                end if
                bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
                bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
              end do
              do m=1,ndim
                bsnew(k)%z(m) = bs(k)%z(m) - dz(m)
              end do

              ! Second child trajectory
              bsnew(nbf+j)%D_big = bs(k)%D_big
              bsnew(nbf+j)%d_pes(in_pes) = (0.0d0,0.0d0)
              do r=1,npes
                if (r.ne.in_pes) then
                  if (x.eq.0) then
                    bsnew(nbf+j)%d_pes(r) = (1.0d0,0.0d0)
                  else
                    bsnew(nbf+j)%d_pes(r) = bs(k)%d_pes(r)
                  end if
                end if
                bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
                bsnew(nbf+j)%a_pes(r) = bsnew(nbf+j)%d_pes(r) * cdexp(i*bsnew(nbf+j)%s_pes(r))
              end do
              do m=1,ndim
                bsnew(nbf+j)%z(m) = bs(k)%z(m) + dz(m)
              end do

            else if (clonetype==3) then

            !!CCS circular distribution of amplitudes - generalised to n potential energy surfaces

            allocate(length(npes),phi(npes))

            !First Child Trajectory

            do r=1,npes
              dummy1 = 0
              do while (dummy1 .eq. 0)
                call random_number(length(r))
                call random_number(phi(r))
                phi(r) = phi(r) * 2.0d0 * pirl
                if (r == 1) then
                  dummy1 = 1
                else
                  do s = 1,r-1
                    !!!! Check that the amplitude values would be sufficiently separated
                    !!!! in the complex plane. This does not seem very efficient however.
                    dist = 0.0d0
                    dist = (length(r)*dcos(phi(r))-length(s)*dcos(phi(s)))**2
                    dist = dist + (length(r)*dsin(phi(r))-length(s)*dsin(phi(s)))**2
                    dist = sqrt(dist)
                    if (dist > qsce) then
                      dummy1 = 1
                    end if
                    !!!! Below checks no two amplitudes would stimulate cloning.
                    !!!! Valid for 2 PESs, but the calculations would change for a more general case
                    !!!! since the cloning is not stimulated by only a single pair of PES components
                    if (dummy1 == 0) then
                      dist = (length(r)**2)*(length(s)**2)/(length(r)**2)+(length(s)**2)
                      if (dist > thresh) then
                        dummy1 = 1
                      end if
                    end if
                  end do
                end if
              end do
            end do

            do r=1,npes
              bsnew(k)%a_pes(r) = cmplx(length(r)*dcos(phi(r)),length(r)*dsin(phi(r)),kind=8)
            end do

            !Second Child Trajectory

            do r=1,npes
              s = npes - (r-1)
              if (r .le. (npes/2.)) then
                bsnew(nbf+j)%a_pes(r) = cmplx(length(s)*dcos(phi(s)),(-1.0d0)*length(s)*dsin(phi(s)),kind=8)
              else
                bsnew(nbf+j)%a_pes(r) = cmplx((-1.0d0)*length(s)*dcos(phi(s)),length(s)*dsin(phi(s)),kind=8)
              end if
              if ((mod(npes,2)==1).and.(r==int(floor(real(npes)/2.))+1)) then
                bsnew(nbf+j)%a_pes(r) = cmplx(0.0d0*dcos(phi(r)),0.0d0*dsin(phi(r)),kind=8)
              end if
            end do

            deallocate (length, phi)

            do m=1,ndim
              bsnew(k)%z(m) = bs(k)%z(m)
              bsnew(nbf+j)%z(m) = bs(k)%z(m)
            end do
            do r= 1,npes
              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
              bsnew(k)%d_pes(r) = bsnew(k)%a_pes(r) * cdexp(-i*bsnew(k)%s_pes(r))
              bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
              bsnew(nbf+j)%d_pes(r) = bsnew(nbf+j)%a_pes(r) * cdexp(-i*bsnew(nbf+j)%s_pes(r))
            end do


  !          !!OAB circular distribution amplitudes

  !            dummy1 = 0
  !            do while (dummy1.eq.0)
  !              dummy2 = 0
  !              do while(dummy2.eq.0)
  !                call random_number(p)
  !                call random_number(q)
  !                pqsqrd = (p**2)+(q**2)
  !                if (pqsqrd <1) then
  !                  dummy2 = 1
  !                end if
  !              end do
  !              p = p/sqrt(pqsqrd)
  !              q = q/sqrt(pqsqrd)
  !              deltaprob = abs((p**2)-(q**2))
  !              if (deltaprob > qsce) then
  !                dummy1 = 1
  !              end if
  !            end do
  !            call random_number(phaseval)
  !            phaseval = phaseval * 2.0d0 * pirl
  !
  !            !Assign First child trajectory a'1
  !            bsnew(k)%a_pes(1) = cmplx(p,0.0d0,kind=8)
  !            !Assign First child trajectory a'2
  !            bsnew(k)%a_pes(2) = cmplx(q*dcos(phaseval),q*dsin(phaseval),kind=8)
  !            !Assign Second child trajectory to be orthogonal a"1
  !            bsnew(nbf+j)%a_pes(1)= cmplx(q*dcos(phaseval),-1.0d0*q*dsin(phaseval),kind=8)
  !            !Assign Second child trajectory to be orthogonal a"2
  !            bsnew(nbf+j)%a_pes(2)= cmplx(-p,0.0d0,kind=8)
  !
  !            do m=1,ndim
  !              bsnew(k)%z(m) = bs(k)%z(m)
  !              bsnew(nbf+j)%z(m) = bs(k)%z(m)
  !            end do
  !            do r= 1,npes
  !              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
  !              bsnew(k)%d_pes(r) = bsnew(k)%a_pes(r) * cdexp(-i*bsnew(k)%s_pes(r))
  !              bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
  !              bsnew(nbf+j)%d_pes(r) = bsnew(nbf+j)%a_pes(r) * cdexp(-i*bsnew(nbf+j)%s_pes(r))
  !            end do

  !            !!CCS square distribution of amplitudes

  !            do r=1,npes
  !              !First child trajectory
  !              sumamps = 0.0d0
  !              q = q * 2.0 - 1.0     ! possibly try modifying the range to enforce close to pes amplitudes
  !              p = p * 2.0 - 1.0
  !              bsnew(k)%d_pes(r) = cmplx(q,p,kind=8)
  !              sumamps = sumamps + dsqrt(dble(cmplx(q,p,kind=8) * cmplx(q,-1.0*p,kind=8)))
  !            end do
  !            do r=1,npes
  !              bsnew(k)%d_pes(r) = bsnew(k)%d_pes(r) / sumamps
  !              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
  !              bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
  !            end do
  !            do m=1,ndim
  !              bsnew(k)%z(m) = bs(k)%z(m)
  !            end do
  !
  !            !Second child trajectory
  !            sumamps = 0.0d0
  !            do r=1,npes
  !              call random_number(q)
  !              call random_number(p)
  !              q = q * 2.0 - 1.0     ! possibly try modifying the range to enforce close to pes amplitudes
  !              p = p * 2.0 - 1.0
  !              bsnew(nbf+j)%d_pes(r) = cmplx(q,p,kind=8)
  !              sumamps = sumamps + dsqrt(dble(cmplx(q,p,kind=8) * cmplx(q,-1.0*p,kind=8)))
  !            end do
  !            do r=1,npes
  !              bsnew(nbf+j)%d_pes(r) = bsnew(nbf+j)%d_pes(r) / sumamps
  !              bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
  !              bsnew(nbf+j)%a_pes(r) = bsnew(nbf+j)%d_pes(r) * cdexp(i*bsnew(nbf+j)%s_pes(r))
  !            end do
  !            do m=1,ndim
  !              bsnew(nbf+j)%z(m) = bs(k)%z(m)
  !            end do

              !Modifications to capital D amplitudes

              bsnew(k)%D_big = bs(k)%D_big * (bsnew(nbf+j)%a_pes(2)*bs(k)%a_pes(1)-bsnew(nbf+j)%a_pes(1)*bs(k)%a_pes(2))
              bsnew(nbf+j)%D_big = bs(k)%D_big * (bsnew(k)%a_pes(1)*bs(k)%a_pes(2)-bsnew(k)%a_pes(2)*bs(k)%a_pes(1))

              bsnew(k)%D_big = bsnew(k)%D_big/&
                      (bsnew(k)%a_pes(1)*bsnew(nbf+j)%a_pes(2)-bsnew(k)%a_pes(2)*bsnew(nbf+j)%a_pes(1))
              bsnew(nbf+j)%D_big = bsnew(nbf+j)%D_big/&
                      (bsnew(k)%a_pes(1)*bsnew(nbf+j)%a_pes(2)-bsnew(k)%a_pes(2)*bsnew(nbf+j)%a_pes(1))

            end if

            write(47756,"(3i5,2es25.17e3)") x, k, nbf+j, abs(bs(k)%a_pes(in_pes)), sqrt(1.-((abs(bs(k)%a_pes(in_pes))**2.0d0)))
            j = j+1

          else

            bsnew(k)%D_big = bs(k)%D_big
            do r=1,npes
              bsnew(k)%d_pes(r) = bs(k)%d_pes(r)
              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
              bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
            end do
            do m=1,ndim
              bsnew(k)%z(m) = bs(k)%z(m)
            end do

          end if

        end do
      else !V1 Cloning moves the basis functions to be cloned 
       
        do k=1,nbf
          bsnew(k)%D_big = bs(k)%D_big
          bsnew2(k)%D_big = bs(k)%D_big
          
          if(clonehere(k)==1) then
            do r=1,npes
              if (r==in_pes)then
                bsnew(k)%d_pes(r)=bs(k)%d_pes(r)
                bsnew2(k)%d_pes(r)=(0.0d0,0.0d0)
              else if (r.ne.in_pes) then
                bsnew2(k)%d_pes(r)=bs(k)%d_pes(r)
                bsnew(k)%d_pes(r)=(0.0d0,0.0d0)
              end if
              bsnew(k)%s_pes(r)= bs(k)%s_pes(r)
              bsnew2(k)%s_pes(r)= bs(k)%s_pes(r)
              bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
              bsnew2(k)%a_pes(r) = bsnew2(k)%d_pes(r) * cdexp(i*bsnew2(k)%s_pes(r))
            end do
          else
            do r=1,npes
              bsnew(k)%d_pes(r) = bs(k)%d_pes(r)
              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
              bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
              bsnew2(k)%d_pes(r) = bs(k)%d_pes(r)
              bsnew2(k)%s_pes(r) = bs(k)%s_pes(r)
              bsnew2(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
            end do
            do m=1,ndim
              bsnew(k)%z(m)=bs(k)%z(m)
              bsnew2(k)%z(m)=bs(k)%z(m)
            end do
          end if
        end do
        call reloc_basis(bsnew2, bs, x)
        call reloc_basis(bsnew, bs, x)
        randomintV1=reps+2000+x
        
        call outbs(bsnew2,randomintV1,mup,muq,time,x)
        write(crrntrep,"(i4.4)") reps
        write(ranomintV1char,"(i4.4)") randomintV1
        call execute_command_line('cp normpop-'//trim(crrntrep)//'.out normpop-'//trim(ranomintV1char)//'.out')
        call deallocbs(bsnew2)

        do k=1,nbfnew
          bs(k)%D_big = bsnew(k)%D_big
          do r=1,npes
            bs(k)%a_pes(r) = bsnew(k)%a_pes(r)
            bs(k)%d_pes(r) = bsnew(k)%d_pes(r)
            bs(k)%s_pes(r) = bsnew(k)%s_pes(r)
          end do
          do m=1,ndim
            bs(k)%z(m) = bsnew(k)%z(m)
          end do
        end do
        call deallocbs(bsnew)
        close (47756)
        return
      end if

      close (47756)

      !Reloc_basis called for original MCEv1 cloning type
      if ((method=="MCEv1")) call reloc_basis(bsnew, bs, x) !.and.(cloneflg.ne."V1"))

      call deallocbs(bs)
      call allocbs(bs, nbfnew)

      do k=1,nbfnew
        bs(k)%D_big = bsnew(k)%D_big
        do r=1,npes
          bs(k)%a_pes(r) = bsnew(k)%a_pes(r)
          bs(k)%d_pes(r) = bsnew(k)%d_pes(r)
          bs(k)%s_pes(r) = bsnew(k)%s_pes(r)
        end do
        do m=1,ndim
          bs(k)%z(m) = bsnew(k)%z(m)
        end do
      end do

      call deallocbs(bsnew)

      n = nbfnew-nbf

      write (6,"(i0,a,i0,a,i0)") sum(clonehere(:)), " bfs cloned in step ", x, &
            ". nbf now = ", nbfnew

      nbf = nbfnew

    end if

    deallocate(clonehere, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating the clone-here array"
      errorflag = 1
      return
    end if

  end subroutine cloning

!***********************************************************************************!
end module bsetalter
