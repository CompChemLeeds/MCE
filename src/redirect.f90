MODULE redirect

  use globvars
  use sb
  use hp
  use mp

!*************************************************************************************************!
!*
!*         Redirection Module
!*           
!*   Contains subroutines for:
!*
!*   1) Redirecting the subroutine which reads in system specific values
!*   2) Redirecting the calculation for the initial z
!*   3) Redirecting the calculation for the single basis function npes x npes Hamiltonian matrix
!*   4) Redirecting the calculation of the derivative of the Hamiltonian
!*      
!*************************************************************************************************!

contains

!--------------------------------------------------------------------------------------------------

  subroutine readparams

    implicit none

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        call readparams_sb
      case ("HP") 
        call readparams_hp
      case ("MP")
        call readparams_mp
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return

  end subroutine readparams    

!--------------------------------------------------------------------------------------------------

  subroutine genzinit(mup, muq, reps)   !   Level 1 Subroutine

    implicit none

    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq
    integer, intent(in) :: reps

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        call genzinit_sb(mup,muq,reps)
      case ("HP")
        call genzinit_hp(mup,muq)
      case ("MP")
        call genzinit_mp(mup,muq)
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return

  end subroutine genzinit

!--------------------------------------------------------------------------------------------------

  subroutine Hord(bs, H, t, reps)
    
    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    type (hamiltonian), dimension (:,:), allocatable, intent(inout) :: H
    real(kind=8), intent (in) :: t
    integer, intent(in) :: reps

    if (errorflag .ne. 0) return

    if (.not.allocated(H)) then
      write(0,"(a)") "Error! Hord has not been allocated in Hord subroutine."    
      errorflag = 1
      return
    end if

    select case (sys)
      case ("SB")
        call Hord_sb(bs, H, t, reps)
      case ("HP")
        call Hord_hp(bs, H, t)
      case ("MP")
        call Hord_mp(bs, H, t)
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select    
    
    return

  end subroutine Hord

!--------------------------------------------------------------------------------------------------

  subroutine Hijdiag(H,z,t, reps)

    implicit none
    complex(kind=8), dimension (:,:), intent(in)::z
    complex(kind=8), dimension(:,:,:), intent (inout)::H
    real(kind=8), intent (in) :: t
    integer, intent(in) :: reps

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        call Hijdiag_sb(H,z,reps)
      case ("HP")
        call Hijdiag_hp(H,z)
      case ("MP")
        call Hijdiag_mp(H,z)
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return   

  end subroutine Hijdiag

!--------------------------------------------------------------------------------------------------

  subroutine dh_dz(dhdz, z, t, reps)

    implicit none
    complex(kind=8),dimension(:,:,:,:), intent(inout) :: dhdz
    complex(kind=8),dimension(:,:),intent(inout)::z 
    real(kind=8), intent (in) :: t 
    integer, intent(in) :: reps

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        dhdz=dh_dz_sb(z,reps)
      case ("HP")
        dhdz=dh_dz_hp(z)
      case ("MP")
        dhdz=dh_dz_mp(z)  
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return

  end subroutine dh_dz

!--------------------------------------------------------------------------------------------------

  subroutine extras(extra, bs)

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    complex(kind=8), intent (inout) :: extra

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        extra=(0.0d0,0.0d0)
      case ("HP")
        extra=disp_hp(bs)
      case ("MP")
        extra=(0.0d0,0.0d0)     
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return

  end subroutine extras

!*************************************************************************************************!

end module redirect
