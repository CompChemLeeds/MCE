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

  subroutine cloning(bs,nbf,x,time,clone, clonenum, reps)
  
    !NOTE: Cloning is only set up for 2 PESs. This needs to be generalised!

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    type(basisfn), dimension(:), allocatable :: bsnew
    real(kind=8), intent(in) :: time
    integer, dimension(:), allocatable, intent(inout) :: clone, clonenum
    integer, intent (inout) :: nbf
    integer, intent (in) :: x, reps
    complex(kind=8), dimension (:), allocatable :: dz
    real(kind=8), dimension (:), allocatable :: dummy_arr
    real(kind=8) :: brforce, normar, sumamps, p, q
    integer, dimension(:), allocatable :: clonehere, clonecopy, clonecopy2 
    integer :: k, m, j, n, nbfnew, ierr, r, clonetype
    character(LEN=3)::rep

    if (errorflag==1) return

    if ((cloneflg=="YES").or.((cloneflg=="BLIND+").and.(x.ne.0))) then
      clonetype = 1 ! 1=conditional cloning, 2=blind cloning
    else if ((cloneflg=="BLIND").or.((cloneflg=="BLIND+").and.(x.eq.0))) then
      clonetype = 2
    else
      write (0,"(2a)") "Cloneflg is invalid. Should be 'YES', 'BLIND' or 'BLIND+', but was ", cloneflg
      errorflag = 1
      return
    end if

    allocate (clonehere(nbf), stat=ierr)
    if (ierr==0) allocate(clonecopy(nbf), stat=ierr)
    if (ierr==0) allocate(clonecopy2(nbf), stat=ierr)
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

    do k=1,nbf
      if (clone(k).lt.x-clonefreq) clone(k)=0
      clonehere(k) = 0
    end do

    ! Build the map of which trajectories to clone, done on the basis function by basis function level (clonetype=1), 
    ! or the entire basis set at once (clonetype=2)

    if (clonetype==1) then
      do k=1,nbf
        normar = 0.0d0
        do r=1,npes
          normar = normar + dconjg(bs(k)%a_pes(r))*bs(k)%a_pes(r)
        end do
        brforce = ((abs(bs(k)%a_pes(1)*bs(k)%a_pes(2))**2.0)/(normar**2.0))
        if ((brforce.gt.thresh).and.(clone(k)==0).and.(clonenum(k).lt.clonemax)) then
          clone(k) = x
          clonehere(k) = 1
        end if 
        clonecopy(k) = clone(k)
        clonecopy2(k) = clonenum(k)
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

    ! Actual cloning section

    if (nbfnew/=nbf) then

      call allocbs(bsnew, nbfnew)

      j=1
      
      write(rep,"(i3.3)") reps
      open(unit=47756,file="Clonetrack-"//trim(rep)//".out",status="old",access="append",iostat=ierr)

      allocate (dz(ndim), stat=ierr)
      if (ierr/=0) then
        write(0,'(a)') "Error in allocating dz array for MCEv1 cloning/relocation"
        errorflag = 1
        return
      end if

!      if (qsc==1) then
!        clonetype = 3
!      else        
      if (method=="MCEv2") then
        clonetype = 1
      else if (method=="MCEv1") then
        clonetype = 2
      end if
      
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
          
            !First child trajectory
            sumamps = 0.0d0
            do r=1,npes
              call random_number(q)
              call random_number(p)
              q = q * 2.0 - 1.0     ! possibly try modifying the range to enforce close to pes amplitudes
              p = p * 2.0 - 1.0
              bsnew(k)%d_pes(r) = cmplx(q,p,kind=8)
              sumamps = sumamps + dsqrt(dble(cmplx(q,p,kind=8) * cmplx(q,-1.0*p,kind=8)))
            end do
            do r=1,npes
              bsnew(k)%d_pes(r) = bsnew(k)%d_pes(r) / sumamps
              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
              bsnew(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
            end do
            do m=1,ndim
              bsnew(k)%z(m) = bs(k)%z(m)
            end do
            
            !Second child trajectory
            sumamps = 0.0d0
            do r=1,npes
              call random_number(q)
              call random_number(p)
              q = q * 2.0 - 1.0     ! possibly try modifying the range to enforce close to pes amplitudes
              p = p * 2.0 - 1.0
              bsnew(nbf+j)%d_pes(r) = cmplx(q,p,kind=8)
              sumamps = sumamps + dsqrt(dble(cmplx(q,p,kind=8) * cmplx(q,-1.0*p,kind=8)))
            end do
            do r=1,npes
              bsnew(nbf+j)%d_pes(r) = bsnew(nbf+j)%d_pes(r) / sumamps
              bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
              bsnew(nbf+j)%a_pes(r) = bsnew(nbf+j)%d_pes(r) * cdexp(i*bsnew(nbf+j)%s_pes(r))
            end do
            do m=1,ndim
              bsnew(nbf+j)%z(m) = bs(k)%z(m)
            end do
            
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
      
      close (47756)
      
      if (method=="MCEv1") call reloc_basis(bsnew, bs, x)
   
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
