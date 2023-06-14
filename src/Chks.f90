MODULE Chks

  use globvars
  use Ham
  use redirect

!***********************************************************************************!
!*
!*         Checking Module
!*           
!*   Contains subroutines for:
!*
!*      1) Checking the initial norm and population sum falls in acceptable range, 
!*         and if desired changes the compression parameter and restarts the basis
!*         set generation.
!*      2) Checking the position components of the coherent state are not too 
!*         widely spaced
!*      3) Checking that the conserved quantites do not deviate too much from their
!*         initial values
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

  subroutine initnormchk(bs,recalcs,restart,alcmprss,absnorm,popsum,trspace)

    implicit none

    type(basisfn), dimension (:), intent (in) :: bs
    complex(kind=8), dimension(:,:), allocatable :: ovrlp
    integer, intent (inout) :: recalcs, restart
    real(kind=8), intent(inout) :: alcmprss, absnorm, popsum
    integer, intent(inout) :: trspace
    complex(kind=8)::normtemp
    integer :: r, k, m, fields, fileun, ierr
    character(len=5) :: loc_basis
    character(LEN=21) :: filenm, filenm2, myfmt
    character(LEN=10):: timestr

    if (errorflag .ne. 0) return

    popsum = 0.0d0

    allocate (ovrlp(size(bs),size(bs)))   
    ovrlp=ovrlpmat(bs)
    
    if (errorflag .ne. 0) return

    normtemp = norm(bs, ovrlp)
    
    if (errorflag .ne. 0) errorflag = 0  ! just in case the initial norm is >1.2

    absnorm = abs(normtemp)

    do r=1,npes
      popsum = popsum + pop(bs, r, ovrlp)
    end do
    deallocate(ovrlp)
    
    if (recalcs.eq.-1) then
      loc_basis = "SWARM"
    else
      loc_basis = basis
    end if
    
    recalcs = 0
 
    if(size(bs).ne.1) then
       
      if (abs(popsum-absnorm).gt.1.0d-10) then
        write(0,"(a)") "Error! Difference between norm and population sum is too high"
        write(0,"(a)") ""
        write(0,"(a,es16.8e3)") "ABS(Norm)  ", absnorm
        write(0,"(a,es16.8e3)") "Popsum     ", popsum
        write(0,"(a,es16.8e3)") "Difference ", abs(popsum-absnorm)
        write(0,"(a)") ""
        restart = 1
      end if

      if (cmprss.eq."N") then      
        if ((absnorm.gt.uplimnorm).or.(absnorm.lt.lowlimnorm)) then
          write(6,"(a,a,e13.5e3)") "Warning. Initial Norm outside established ",&
                       "parameters, with a value of ", absnorm
          write(6,"(a)") ""
          if ((loc_basis.eq."SWARM").or.(loc_basis.eq."SWTRN")) restart = 1
        end if
      else
        if (absnorm.lt.lowlimnorm) then
          write(6,'(a,es16.8e3)') "Initial Norm too low with a value of ", absnorm
          if (loc_basis.eq."SWTRN") then
            write(6,"(a,i0)") "Reducing train spacing spacing to ", int(real(trspace) * 0.95)
            write(6,"(a)") ""
            trspace = int(real(trspace) * 0.95)
          else if (loc_basis.eq."SWARM") then
            write(6,"(a,es16.8e3)") "Increasing compression parameter to", 1/(alcmprss * 0.95d0)
            write(6,"(a)") ""
            alcmprss = alcmprss * 0.95d0
          end if 
          restart = 1
        else if (absnorm.gt.uplimnorm) then
          write(6,'(a,es16.8e3)') "Initial Norm too high with a value of ", absnorm
          if (loc_basis.eq."GRID") then
            write(6,"(a,i0)") "Increasing train spacing to ", int(real(trspace) * 1.05)
            write(6,"(a)") ""
            trspace = int(real(trspace) * 1.05)
          else if (loc_basis.eq."SWARM") then
            write(6,"(a,es16.8e3)") "Reducing compression parameter to", 1/(alcmprss * 1.05d0)
            write(6,"(a)") ""
            alcmprss = alcmprss * 1.05d0
          end if 
          restart = 1
        end if
      end if
      
      if ((restart.eq.1).and.(recalcs.lt.Ntries)) then
        if (cmprss.eq."N") then
          recalcs = recalcs + 1
        end if
        write(6,"(a)") "Recalculating..."
        write(6,"(a)") ""
        
        write(timestr,"(i3.3)") recalcs
        
        return
      else
        return
      end if
   
    else
   
      return
   
    end if

  end subroutine initnormchk

!------------------------------------------------------------------------------------

  subroutine enchk(bf,t,n,redo,k,reps)

    implicit none
    type(basisfn), intent(inout) :: bf
    real(kind=8), intent(in) :: t
    integer, intent(inout) :: n, redo
    integer, intent(in) :: k, reps
    
    complex(kind=8),dimension(:,:,:), allocatable::H
    complex(kind=8),dimension(:,:), allocatable::z
    real(kind=8)::Echk
    integer :: j,r,s,ierr

    if (errorflag/=0) return

    if (ECheck.eq."YES") then
    
      allocate(H(1,npes,npes), stat = ierr)
      if (ierr==0) allocate (z(1,ndim), stat = ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error in H allocation in genbasis"
        errorflag=1
        return
      end if
      
      call Hijdiag(H, z, t, reps)
      
      Echk = 0.0d0
      
      do j=1,size(H,1)
        do r=1,size(H,2)
          do s=1,size(H,3)
            Echk = Echk + H(j,r,s)
          end do
        end do
      end do
      
      if ((Echk.gt.Ebfmax).or.(Echk.lt.Ebfmin)) then
        if (n.lt.Ntries) then
          n = n+1
          write(6,"(a)")"Basis ", k, " did not meet energy requirements. ",&
                 "Recalculating..."
          redo=1
        else
          write(6,"(a)")"Basis ", k, " recalculated ", n, "times but still outside ",&
                  "acceptable region."
          write(0,"(a)") "Multiple recalculations did not find an adequate basis"
          write(0,"(a)") "Terminating calculation"
          errorflag = 1
          redo=0
        end if
      else
        redo=0
      end if
    else
      redo=0
    end if

    if ((redo/=1).and.(redo/=0)) then
      write(0,"(a)") "Error! Somehow, the redo flag is not 1 or 0"
      errorflag = 1
      return
    end if

    return

  end subroutine enchk

!------------------------------------------------------------------------------------

  subroutine trajchk(bs)   !   Level 1 Subroutine

    implicit none

    type(basisfn),  dimension(:), intent(in) :: bs
    complex(kind=8), dimension(ndim) :: z
    complex(kind=8) :: trajq
    integer :: j, m, flag

    if (errorflag .ne. 0) return

    do j = 1,size(bs)
      flag = 0
      z = bs(j)%z
      do m=1,ndim
        trajq = dble(z(m))*dsqrt(2.0d0/gam)
        if (abs(trajq).gt.20000) then
          write(6,'(a,i0,a,i0,a,e16.8)') "Trajectory ", j, " in dof ", m, &
                        " is equal to ", dble(z(m))*dsqrt(2.0d0)
          flag = flag + 1
        end if
      end do
    end do

    if (flag.gt.0) then
      errorflag = 1
      write(6,"(i0,a,a)") flag, " trajectories have position components outside ",&
                      "acceptable range."
      return
    end if

    return

  end subroutine trajchk

!------------------------------------------------------------------------------------

  subroutine conservchk(initehr, initnorm, absehr, absnorm, reps)   !   Level 1 Subroutine

    implicit none

    real(kind=8), intent (in) :: initehr, initnorm, absehr, absnorm
    character(LEN=15) :: filenm
    integer, intent(in) :: reps
    character(LEN=3):: rep

    if (errorflag .ne. 0) return

    write(rep,"(i3.3)") reps

    filenm = "normpop-"//trim(rep)//".out"

    if (abs(initnorm-absnorm).ge.1.0d-3) then
      write (0,"(a)") ""
      write (0,"(a)") "*************************Simulation Failed*************************"
      write (0,"(a)") "*************Norm deviated too much from initial value*************"
      errorflag = 1
      return
    end if

    if ((abs(1.0d0-(absehr/initehr)).ge.1.0d-2).and.(method=="MCEv2")) then
      write (0,"(a)") ""
      write (0,"(a)") "*************************Simulation Failed*************************"
      write (0,"(a)") "*******Ehrenfest Energy deviated too much from initial value*******"
      errorflag = 1
      return
    end if

    return

  end subroutine conservchk    

!*************************************************************************************************!  

end module chks


