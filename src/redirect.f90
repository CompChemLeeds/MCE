MODULE redirect

  use globvars
  use sb
  use vp
  use hp
  use fp
  use mp
  use iv
  use cp
  use hh

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
      case ("VP")
        call readparams_vp
      case ("HP") 
        call readparams_hp
      case ("FP")
        call readparams_fp
      case ("MP")
        call readparams_mp
      case ("IV")
        call readparams_iv
      case ("CP")
        call readparams_cp
      case ("HH")
        call readparams_hh
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return

  end subroutine readparams    

!--------------------------------------------------------------------------------------------------

  subroutine genzinit(mup, muq)   !   Level 1 Subroutine

    implicit none

    real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        call genzinit_sb(mup,muq)
      case ("VP")
        call genzinit_vp(mup,muq)
      case ("HP")
        call genzinit_hp(mup,muq)
      case ("FP")
        call genzinit_fp(mup,muq)
      case ("MP")
        call genzinit_mp(mup,muq)
      case ("IV")
        call genzinit_iv(mup,muq)
      case ("CP")
        call genzinit_cp(mup,muq)
      case ("HH")
        call genzinit_hh(mup,muq)
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return

  end subroutine genzinit

!--------------------------------------------------------------------------------------------------

  subroutine Hord(bs, H, t)
    
    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    type (hamiltonian), dimension (:,:), allocatable, intent(inout) :: H
    real(kind=8), intent (in) :: t

    if (errorflag .ne. 0) return

    if (.not.allocated(H)) then
      write(0,"(a)") "Error! Hord has not been allocated in Hord subroutine."    
      errorflag = 1
      return
    end if

    select case (sys)
      case ("SB")
        call Hord_sb(bs, H, t)
      case ("VP")
        call Hord_vp(bs, H, t)
      case ("HP")
        call Hord_hp(bs, H, t)
      case ("FP")
        call Hord_fp(bs, H, t)
      case ("MP")
        call Hord_mp(bs, H, t)
      case ("IV")
        call Hord_iv(bs, H, t)
      case ("CP")
        call Hord_cp(bs, H, t)
      case ("HH")
        call Hord_hh(bs, H, t)
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select    
    
    return

  end subroutine Hord

!--------------------------------------------------------------------------------------------------

  subroutine Hijdiag(H,z,t)

    implicit none
    complex(kind=8), dimension (:,:), intent(in)::z
    complex(kind=8), dimension(:,:,:), intent (inout)::H
    real(kind=8), intent (in) :: t

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        call Hijdiag_sb(H,z)
      case ("VP")
        call Hijdiag_vp(H,z)
      case ("HP")
        call Hijdiag_hp(H,z)
      case ("FP")
        call Hijdiag_fp(H,z)
      case ("MP")
        call Hijdiag_mp(H,z)
      case ("IV")
        call Hijdiag_iv(H,z,t)
      case ("CP")
        call Hijdiag_cp(H,z,t)
      case ("HH")
        call Hijdiag_hh(H,z)
      case default
        write(0,"(a)") "Error! The system was not recognised!"
        write(0,"(a)") "If you are seeing this something is terribly wrong"
        errorflag=1
    end select

    return   

  end subroutine Hijdiag

!--------------------------------------------------------------------------------------------------

  subroutine dh_dz(dhdz, z, t)

    implicit none
    complex(kind=8),dimension(:,:,:,:), intent(inout) :: dhdz
    complex(kind=8),dimension(:,:),intent(inout)::z 
    real(kind=8), intent (in) :: t 

    if (errorflag .ne. 0) return

    select case (sys)
      case ("SB")
        dhdz=dh_dz_sb(z)
      case ("VP")
        dhdz=dh_dz_vp(z)
      case ("HP")
        dhdz=dh_dz_hp(z)
      case ("FP")
        dhdz=dh_dz_fp(z)
      case ("MP")
        dhdz=dh_dz_mp(z)  
      case ("IV")
        dhdz=dh_dz_iv(z,t)
      case ("CP")
        dhdz=dh_dz_cp(z,t)    
      case ("HH")
        dhdz=dh_dz_hh(z)   
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
      case ("VP")
        extra=disp_vp(bs)
      case ("HP")
        extra=disp_hp(bs)
      case ("FP")
        extra=(0.0d0,0.0d0)
      case ("MP")
        extra=(0.0d0,0.0d0)  
      case ("IV")
        extra=dipole_iv(bs)
      case ("CP")
        extra=dipole_cp(bs)
      case ("HH")
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
