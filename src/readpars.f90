MODULE readpars

  use globvars
  use alarrays
  use Ham
  use redirect

!*************************************************************************************************!
!*
!*         Input File Reading Module
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the run conditions (gen,prop,cmprss,method,reptot,conjflg)
!*      2) Reading the system type (currently only Spin Boson supported)
!*      3) Reading the energy cutoff parameters (ECheck,Ntries,Ebfmin,Ebfmax)
!*      4) Reading the basis set parameters (ndim,in_nbf,matfun,npes,in_pes,grid)
!*      5) Reading initial wavefunction parameters (initialcmprss,gam,mu,hbar and calculation of sigp & sigq)
!*      6) Reading a pre-calculated basis set (inc. bs params and all bs values)
!*      7) Reading time propagation parameters (dtmin,dtmax,dtinit,timeend,timestrt,step)
!*      
!*************************************************************************************************!

contains

!*************************************************************************************************!
!          Shared Reading Subroutines
!*************************************************************************************************!

  subroutine readrunconds   !   Level 1 Subroutine

    implicit none
    character(LEN=100)::LINE, LINE2
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=140, file='input.dat', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening runconds.dat file'
      errorflag = 1
      return
    end if

    read(140,*,iostat=ierr)LINE

    do while (ierr==0)

      if (LINE=='gen') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,LINE2
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading basis set generation status"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          gen = "Y"
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          gen = "N"
        else
          write(0,"(a,a)") "Error. gen value must be YES/NO. Read ", trim(LINE2)
        end if
        n=n+1
      else if (LINE=='prop') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,LINE2
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading basis set propagation status"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          prop = "Y"
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          prop = "N"
        else
          write(0,"(a,a)") "Error. prop value must be YES/NO. Read ", trim(LINE2)
        end if
        n=n+1   
      else if (LINE=='restart') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,restrtflg
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading basis set propagation status"
          errorflag = 1
          return
        end if
        if ((restrtflg.ne.0).and.(restrtflg.ne.1)) then
          write(0,"(a,a)") "Error. Restart flag must be 1/0. Read ", trim(LINE2)
        end if
        n=n+1             
      else if (LINE=='cmprss') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,LINE2
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading compression parameter change status"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          cmprss = "Y"
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          cmprss = "N"
        else
          write(0,"(a,a)") "Error. cmprss value must be YES/NO. Read ", trim(LINE2)
        end if
        n=n+1 
      else if (LINE=='method') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,LINE2
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading basis set propagation method"
          errorflag = 1
          return
        end if
        if ((LINE2(1:5).eq.'mcev1').or.(LINE2(1:5).eq.'MCEv1')) then
          method = "MCEv1"
        else if ((LINE2(1:5).eq.'mcev2').or.(LINE2(1:5).eq.'MCEv2')) then
          method = "MCEv2"
        else if ((LINE2(1:3).eq.'ccs').or.(LINE2(1:3).eq.'CCS')) then
          method = "CCS"
        else
          write(0,"(a,a,a)") "Error. Method must be MCEv1, MCEv2, or CCS. Read ", trim(LINE2)
        end if
        n=n+1
      else if (LINE=='Repeats') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,reptot
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading number of repeats"
          errorflag = 1
          return
        end if
        n=n+1                    
      else if (LINE=='Conjugate_Repeats') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,LINE2
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading conjugate repeats flag"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          conjflg = 1
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          conjflg = 0
        else
          write(0,"(a,a)") "Error. Conjugate repeats flag must be Yes or No. Read ", trim(LINE2)
        end if
        n=n+1                     
      end if

      read(140,*,iostat=ierr) LINE

    end do

    close(140)

    if ((gen.eq."N").and.(prop.eq."N")) then
      write(0,"(a)") "Error! Run conditions are for no basis set generation or propagation. So what now genius?"
      errorflag=1
      return
    end if

    if ((conjflg==1).and.(mod(reptot,2).ne.0)) then
      write(6,"(a)")"Warning! An odd number of repeats selected but conjugate repetition chosen!"
      write(6,"(a)")"Incrementing repeat total to even number"
      reptot = reptot + 1
    end if

    if ((conjflg==1).and.(gen=="N")) then
      write(0,"(a)")"Error! Conjugate repetition is not compatible for simulations with pre-calculated basis set"
      errorflag=1
      return
    end if
    
    if ((conjflg==1).and.(restrtflg==1)) then
      write(0,"(a)")"Error! Conjugate repetition is not compatible for restarted simulations"
      errorflag=1
      return
    end if

    if (n.ne.7) then
      write(0,"(a)") "Not all required variables read in readrunconds subroutine"
      write(0,"(a,i0,a)") "Read a total of ", n, "of an expected 8 parameters"
      errorflag = 1
      return
    end if

    return

  end subroutine readrunconds

!--------------------------------------------------------------------------------------------------

  subroutine readsys   !   Level 1 Subroutine

    implicit none
    character(LEN=100)::LINE
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=127, file='input.dat', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening input.dat file'
      errorflag = 1
      return
    end if

    read(127,*,iostat=ierr)LINE

    do while (ierr==0)

      if(LINE=='System:') then
        backspace(127)
        read(127,*,iostat=ierr)LINE,sys
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading System Name"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='SpecDen') then
        backspace(127)
        read(127,*,iostat=ierr)LINE,specden
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Spin Boson Spectral Density name"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='freqflg') then
        backspace(127)
        read(127,*,iostat=ierr)LINE,freqflg_sb
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Spin Boson Frequency calculation flag"
          errorflag = 1
          return
        end if
        n = n+1
      end if
      read(127,*,iostat=ierr) LINE

    end do

    close(127)
    
    if ((sys.ne."SB").and.(sys.ne."HP").and.(sys.ne."MP")) then
      write(0,"(2a)") "System not recognised. Please recheck input.dat file. Read value of ", sys
      errorflag = 1
      return
    end if
    
    if ((sys.eq."SB").and.(specden.ne."EXP").and.(specden.ne."DL").and.(specden.ne."UBO").and.(specden.ne."LHC")) then
      write(0,"(2a)") "Spectral Density not recognised. Must be EXP, DL, UBO or LHC, but read ", specden
      errorflag = 1
      return
    end if
    
    if (specden.eq."LHC") then
      write(0,"(a)") "The LHC-II Spectral Density calculation is not full set up yet. Please change the spectral density"
      errorflag = 1
      return
    end if
      
    if ((freqflg_sb.eq.0).and.((specden=="UBO").or.(specden=="LHC"))) then 
      write(0,"(3a)") "Frequency flag set for self calculation. This is not compatible with the ", specden, &
                     & " spectral density, which only works with pre-calculation. Please revise the running parameters"
      errorflag = 1
      return
    end if

    if (n.ne.3) then
      write(0,"(a)") "Not all required variables read in readsys subroutine"
      errorflag = 1
      return
    end if  

    call readparams

    return

  end subroutine readsys

!--------------------------------------------------------------------------------------------------

  subroutine readecut   !   Level 1 Subroutine

    implicit none
    character(LEN=100)::LINE, LINE2
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    open(unit=128, file='inham.dat', status='old', iostat=ierr)

    if (ierr.ne.0) then
      write(0,"(a)") 'Error in opening inham.dat file'
      errorflag = 1
      return
    end if

    read(128,*,iostat=ierr)LINE

    do while (ierr==0)

      if(LINE=='ECheck') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,LINE2
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading ECheck value"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          Echeck = "YES"
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          Echeck = "NO"
        else
          write(0,"(a,a)") "Error. Echeck value must be YES/NO. Read ", trim(LINE2)
        end if
        n = n+1
      else if(LINE=='Ntries') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,Ntries
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Ntries value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='Ebfmin') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,Ebfmin
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Ebfmin value"
          errorflag = 1
          return
        end if
        n = n+1
      else if(LINE=='Ebfmax') then
        backspace(128)
        read(128,*,iostat=ierr)LINE,Ebfmax
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading Ebfmax value"
          errorflag = 1
          return
        end if
        n = n+1
      end if

      read(128,*,iostat=ierr) LINE

    end do

    close (128)
    
    if (Ebfmin.ge.Ebfmax) then
      write(0,"(a)") "Invalid values for Ebfmin and/or Ebfmax. Max must be higher than min"
      errorflag = 1
      return
    end if

    if (n.ne.4) then
      write(0,"(a)") "Not all required variables read in readecut subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine readecut


!*************************************************************************************************!
!          Reading Subroutines for Basis Set Generation
!*************************************************************************************************!

  subroutine readbsparams   !   Level 1 Subroutine

    IMPLICIT NONE
    character(LEN=100)::LINE, LINE2
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    OPEN(UNIT=127, FILE='input.dat',STATUS='OLD', iostat=ierr)

    if (ierr .ne. 0) then
      write(0,"(a)") 'Error in opening input.dat file'
      errorflag = 1
      return
    end if

    read(127,*,iostat=ierr)LINE

    do while (ierr==0)

      if(LINE== "ndim") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,ndim 
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading ndim"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="in_nbf") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,in_nbf
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading in_nbf"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="matfun") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,matfun
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading matrix function"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="randfunc") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,randfunc
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading random number function"
          call flush(0)
          errorflag = 1
          return
        end if
        call flush(6)
        n = n+1
      else if (LINE=="npes") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,npes
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading npes"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="in_PES") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,in_pes
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading in_PES"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="basis") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,basis
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading basis option"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="trainsp") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,trainsp
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading train spacing"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="train_len") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,train_len
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading default number of basis functions per train."
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="swtrn_swarm_size") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,swtrn_swrm
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading size of central SWTRN swarm."
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="Cloning") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,LINE2
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading cloning flag"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          cloneflg = "YES"
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          cloneflg = "NO"
        else if ((LINE2(1:1).eq.'b').or.(LINE2(1:1).eq.'B')) then
          if (LINE2(6:6).eq.'+') then
            cloneflg = "BLIND+"
          else
            cloneflg = "BLIND"
          end if
        else
          write(0,"(a,a)") "Error. cloneflg value must be YES/NO/BLIND/BLIND+. Read ", trim(LINE2)
        end if
        n = n+1
      else if (LINE=="Threshold") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,thresh
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading cloning threshold"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="max_cloning") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,clonemax
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading maximum number of clones"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="clon_freq") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,clonefreq
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading cloning frequency"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="qss") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,qss
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading qss amplitudes"
          errorflag = 1
          return
        end if
        n = n+1
      end if

      read(127,*,iostat=ierr) LINE

    end do

    close(127) 
    
    if ((in_pes.gt.npes).or.(in_pes.le.0)) then
      write(0,"(a)") "Initial PES does not exist"
      ierr=-1
      errorflag = 1
      return
    end if

    if (in_nbf.le.0) then
      write(0,"(a)") "Number of Basis Functions <= 0"
      ierr=-1
      errorflag = 1
      return
    end if

    if (npes.le.0) then
      write(0,"(a)") "Number of PES <= 0"
      ierr=-1
      errorflag = 1
      return
    end if

    if ((npes.eq.1).and.(trim(method)/="CCS")) then
      write(0,"(a)") "Only one PES selected, but propagation method is not CCS!"
      ierr=-1
      errorflag=1
      return
    end if

    if (ndim.le.0) then
      write(0,"(a)") "Number of Degrees of Freedom <= 0"
      ierr=-1
      errorflag = 1
      return
    end if
    
    if (qss==1) then
      if (npes.ne.2) then
        write(0,"(a)") "qss amplitudes are only valid with npes = 2."
        write(0,"(a,i0)") "Npes was read as ", npes
        ierr=-1
        errorflag=1
        return
      end if
    else if (qss.ne.0) then
      write(0,"(a)") "Only allowed values for the qss flag are 1 or 0"
      write(0,"(a,i0)") "Value read was ", qss
      ierr=-1
      errorflag = 1
      return
    end if    
    
    if ((basis.ne.'SWARM').and.(basis.ne.'swarm').and.(basis.ne.'SWTRN').and.(basis.ne.'swtrn')) then
      write(0,"(a,a)") "Invalid value for basis. Must be SWARM/SWTRN and all upper/lower case. Value is ", basis
      errorflag = 1
      return
    else if ((basis.eq.'SWARM').or.(basis.eq.'swarm'))then
      basis = 'SWARM'
    else if ((basis.eq.'SWTRN').or.(basis.eq.'swtrn'))then
      basis = 'SWTRN'
    end if
    
    if (thresh .lt. 0.05d0) then
      write(0,"(a)") "Cloning threshold too small. Setting to default value of 0.249"
      thresh = 0.249d0
    else if (thresh .ge. 0.25d0) then
      write(0,"(a)") "Cloning threshold larger than allowed limits. Setting to default value of 0.249"
      thresh = 0.249d0
    end if

    if ((matfun.ne.'zgesv').and.(matfun.ne.'ZGESV').and.(matfun.ne.'zheev').and.(matfun.ne.'ZHEEV')) then
      write(0,"(a,a)") "Invalid value for matrix function. Must be ZGESV/zgesv or ZHEEV/zheev. Value is ", matfun
      errorflag = 1
      return
    else if ((matfun.eq.'zgesv').or.(matfun.eq.'ZGESV')) then
      matfun = 'zgesv'
    else
      matfun = 'zheev'
    end if
    
    if ((randfunc.ne.'ZBQL').and.(randfunc.ne.'zbql').and.(randfunc.ne.'gaus').and.(randfunc.ne.'GAUS')) then
      write(0,"(a,a)") "Invalid value for random number function. Must be ZBQL/zbql or GAUS/gaus. Value is ", randfunc
      errorflag = 1
      return
    else if ((randfunc.ne.'ZBQL').or.(randfunc.ne.'zbql')) then
      randfunc = 'ZBQL'
    else
      randfunc = 'GAUS'
    end if      
    
    if (n.ne.15) then
      write(0,"(a,i0)") "Not all required variables read in readbsparams subroutine. n=", n
      errorflag = 1
      return
    end if

    return

  end subroutine readbsparams
  
!--------------------------------------------------------------------------------------------------

  subroutine checkparams    !   Subroutiner checks that parameters are compatible with each other
  
    implicit none
    
    if (errorflag.ne.0) return
    
    !!!!!! Basis specific checks !!!!!!!
    
    if (basis.eq."SWTRN") then
      if (step.eq."A") then
        write(0,"(a)") "Trains are only valid for static stepsizes."
        errorflag = 1
        return
      end if
      if (trainsp.le.0) then
        write(0,"(a)") "Error! Spacing between basis set train elements is <=0"
        errorflag = 1
        return
      end if
      in_nbf = swtrn_swrm * train_len
      write(6,"(a,i0)") "Setting the total initial size of the basis set to ", in_nbf
      if ((method.ne."MCEv2").and.(method.ne."CCS")) then
        write(0,"(a)") "Error! Trains can only work with MCEv2 or CCS."
        write(0,"(a)") "MCEv1 cannot calculate the amplitudes correctly"
        errorflag = 1
        return
      end if        
    end if
       
    !!!!!! Basis set change parameter check!!!!!!

    if (cloneflg.ne.'NO') then
      if ((method.ne."MCEv2").and.(method.ne."MCEv1")) then
        write (0,"(a)") "Cloning can only work with MCEv2 or MCEv1"
        errorflag = 1
        return
      else if ((method.eq."MCEv2").and.(basis.ne."SWARM").and.(basis.ne."SWTRN")) then
        write (0,"(a)") "Cloning with MCEv2 can only work on a swarm or a swarm of trains"
        errorflag = 1
        return
      else if ((method.eq."MCEv2").and.(basis.ne."SWARM")) then
        write (0,"(a)") "Cloning with MCEv1 can only work on a swarm"
        errorflag = 1
        return        
      end if
      if (clonemax.lt.1) then
        write(0,"(a)") "Maximum number of clones allowed is less than 1 but cloning is enabled!"
        write(0,"(a)") "These conditions are incompatible"
        errorflag = 1
        return
      else if (clonemax.gt.20) then
        write(0,"(a)") "Maximum number of clones allowed is 20!"
!        write(0,"(a)") "This will mean that over 1024 clones could be made from EACH member of the initial basis set"
        write(0,"(a)") "Try again with a lower number (8 should be enough for anyone)"
        errorflag = 1
        return
      end if
      if (sys.ne."SB") then
        write(0,"(a)") "Cloning can only be performed on multi-PES systems, currently SB, or VP"
        errorflag = 1
        return
      end if
      if (clonefreq.lt.50) then
        write(0,"(a)") "A Cloning frequency of less than 50 could result in too rapid cloning at the beginning of the simulation"
        errorflag = 1
        return
      else if (clonefreq.gt.1000) then
        write(0,"(a)") "A cloning frequency of greater than 1000 will result in very little effect of the cloning procedure"
        errorflag = 1
        return
      end if
    end if 
        
   
    !!!!!!!! Check the System!!!!!!  
    
    select case (sys)
      case ("SB")
        if (npes.lt.2) then
          write(0,"(a)") "Spin Boson model must have at least 2 pes'"
          errorflag = 1
          return
        end if
        if ((method.ne."MCEv1").and.(method.ne."MCEv2")) then
          write(0,"(a)") "Spin Boson model can only be simulated by MCEv1 or MCEv2"
          errorflag = 1
          return
        end if
        if ((basis.ne."SWARM").and.(basis.ne."SWTRN")) then
          write(0,"(a)") "This method must be used with Swarms or Swarm-Trains."
          errorflag = 1
          return
        end if 
      case ("HP")
        if (npes.ne.1) then
          write(0,"(a)") "Harmonic Potential only valid for 1 PES"
          errorflag = 1
          return
        end if
        if (trim(method).ne."CCS") then
          write(0,"(a)") "Harmonic Potential can only be simulated by CCS"
          errorflag = 1
          return  
        end if
      case ("MP")
        if (npes.ne.1) then
          write(0,"(a)") "Morse Potential only valid for 1 PES"
          errorflag = 1
          return
        end if
        if (trim(method).ne."CCS") then
          write(0,"(a)") "Morse Potential can only be simulated by CCS"
          errorflag = 1
          return  
        end if 
      case default
        write(0,"(a)") "System is not recognised. Value is ", sys
        errorflag = 1
        return
    end select  
    
  end subroutine checkparams

!--------------------------------------------------------------------------------------------------

  subroutine readzparams   !   Level 1 Subroutine

    IMPLICIT NONE
    character(LEN=100)::LINE
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    OPEN(UNIT=127, FILE='input.dat',STATUS='OLD', iostat=ierr)

    if (ierr .ne. 0) then
      write(0,"(a)") 'Error in opening input.dat file'
      errorflag = 1
      return
    end if

    read(127,*,iostat=ierr)LINE

    do while (ierr==0)

      if (LINE=="ALCMP") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,initalcmprss
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading compression parameter"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="gamma") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,gam
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading gamma factor"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="mu") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,mu
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading mu (centre of initial random gaussian)"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="hbar") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,hbar
        if(ierr.ne.0) then
          write(6,"(a)")  "Error reading hbar value. This value is optional, defaulting to 1"
          hbar = 1.0d0
        else if (hbar.ne.1.0d0) then 
          write(6,"(a,es12.5)") "hbar changed from default to ", hbar
        end if
        n = n+1
      end if

      read(127,*,iostat=ierr) LINE

    end do

    close(127)

    sigp = sqrt(1.0d0/(2.0d0*gam))
    sigq = sqrt(gam/2.0d0)

    if (n.ne.4) then
      write(0,"(a,i0)") "Not all required variables read in readzparams subroutine. n = ", n
      errorflag = 1
      return
    end if

    return

  end subroutine readzparams

!*************************************************************************************************!
!          Reading Subroutines for Basis Set Propagation
!*************************************************************************************************!

  subroutine readbasis(bs, mup, muq, rep, t, nbf)   !   Level 1 Subroutine

    implicit none
    type(basisfn), dimension (:), allocatable, intent(inout) :: bs
    integer::ierr, n, j, k, m, r, cflg, bsunit
    real(kind=8), dimension(:), allocatable, intent(out) :: mup, muq 
    real(kind=8), intent(inout) :: t
    integer, intent(inout) :: nbf
    integer, intent(in) :: rep  
    character(LEN=100)::LINE
    character(LEN=14)::filename
    real(kind=8)::rl, im
    complex(kind = 8) :: dsum1

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0
    cflg = 0

    write(6,"(a)")"Starting read subroutine"
    call flush(6)

    write(filename,"(a,i4.4,a)") "Outbs-", rep, ".out"

    write(6,"(a,a)") "Opening file ", trim(filename)
    call flush(6)
    
    bsunit = 200+rep

    open(unit=bsunit, file=filename, status="old", iostat=ierr)

    if (ierr .ne. 0) then
      write(0,"(3a)") 'Error in opening ', trim(filename),' file'
      call flush(0)
      errorflag = 1
      return
    end if

    read(bsunit,*,iostat=ierr)LINE
    call flush(6)

    do while ((LINE.ne."zinit").and.(ierr==0))
      if (LINE=="ndof") then
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,ndim
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading ndim"
          call flush(0)
          errorflag = 1
          return
        end if
        write(6,"(a,i0)") "ndim    = ", ndim
        call flush(6)
        n = n+1
      else if (LINE=="nconf") then
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,npes
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading npes"
          call flush(0)
          errorflag = 1
          return
        end if
        write(6,"(a,i0)") "npes    = ", npes
        call flush(6)
        n = n+1
      else if (LINE=="nbasisfns") then
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,nbf
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading nbf"
          call flush(0)
          errorflag = 1
          return
        end if
        write(6,"(a,i0)") "nbf     = ", nbf
        call flush(6)
        n = n+1
      else if (LINE=="initial_PES") then
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,in_pes
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading in_PES"
          call flush(0)
          errorflag = 1
          return
        end if
        write(6,"(a,i0)") "in_pes  = ", in_pes
        call flush(6)
        n = n+1
      else if (LINE=="matfun") then
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,matfun
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading matrix function"
          call flush(0)
          errorflag = 1
          return
        end if
        write(6,"(a,a)") "matfun  = ", matfun
        call flush(6)
        n = n+1
      else if (LINE=="time") then
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,t
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading time"
          call flush(0)
          errorflag = 1
          return
        end if
        write(6,"(a,es16.8e3)") "time =", t
        call flush(6)
        n = n+1
      end if
      read(bsunit,*,iostat=ierr)LINE
    end do
    
    call flush(6)
    call flush(0)

    if (n.ne.6) then
      write(0,"(a,i0,a)") "Error in reading parameters. Only ", n, " of 6 parameters read."
      errorflag = 1
      return
    end if 
    
    call flush(6)
    call flush(0)

    allocate (mup(ndim), stat=ierr)
    if (ierr == 0) allocate (muq(ndim), stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error in allocation of mup and muq"
      errorflag=1
    end if

    if (LINE.ne."zinit") then
      write(0,"(a,a)") "Error! Expected zinit, but read ", trim(LINE)
    end if
    backspace(bsunit)

    do m=1,ndim
      read(bsunit,*,iostat=ierr)LINE, j, muq(j), mup(j)
      if(ierr.ne.0) then
        write(0,"(a,a)")  "Error reading zinit value ", m
        errorflag = 1
        return
      end if
      if(m.ne.j) then
        write(0,"(a,a)")  "Error! Count mismatch in zinit. Expected ", m, "but read ", j
        errorflag = 1
        return
      end if
    end do
     
    read(bsunit,*,iostat=ierr)LINE

    if (LINE.ne."basis") then
      write(0,"(a,a)") "Error! Expected basis, but read ", trim(LINE)
    end if

    backspace(bsunit)

    if (size(bs).ne.nbf) then
      write(6,"(a)") "Basis set size has changed. Reallocating basis set."
      call deallocbs(bs)
      call allocbs(bs, nbf)
    end if

    do j=1,nbf
      read(bsunit,*,iostat=ierr)LINE,k
      if(k.ne.j) then
        write(0,"(a,i2,a,i2)") "Error. Expected basis function ", j, " but got ", k
      end if
      read (bsunit,*,iostat=ierr)LINE
      if (LINE.ne."D") then
        write(0,"(a,a)") "Error! Expected D but read ", trim(LINE)
      end if
      backspace(bsunit)
      read(bsunit,*,iostat=ierr)LINE,rl, im
      if (LINE.eq."D") then
        bs(j)%D_big=cmplx(rl,im,kind=8)
      else
        write(0,"(a)") "Something has gone very wrong here"
        errorflag = 1
        return
      end if
      do r=1,npes
        read(bsunit,*,iostat=ierr)LINE
        if (LINE.ne."a") then
          write(0,"(a,a)") "Error! Expected a but read ", trim(LINE)
        end if
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,k,rl, im
        if (k.ne.r) then
          write (0,"(a,i2,a,i2)") "Error. Expected a from pes ", r, "but got ", k
        end if
        bs(j)%a_pes(r)=cmplx(rl,im,kind=8)
      end do
      do r=1,npes
        read(bsunit,*,iostat=ierr)LINE
        if (LINE.ne."d") then
          write(0,"(a,a)") "Error! Expected d but read ", trim(LINE)
        end if
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,k,rl, im
        if (k.ne.r) then
          write(0,"(a,i2,a,i2)") "Error. Expected d from pes ", r, "but got ", k
        end if
        bs(j)%d_pes(r)=cmplx(rl,im,kind=8)
      end do
      do r=1,npes
        read(bsunit,*,iostat=ierr)LINE
        if (LINE.ne."s") then
          write(0,"(a,a)") "Error! Expected s but read ", trim(LINE)
        end if
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,k,rl
        if (k.ne.r) then
          write(0,"(a,i2,a,i2)") "Error. Expected s from pes ", r, "but got ", k
        end if
        bs(j)%s_pes(r)=rl
      end do
      do m=1,ndim
        read(bsunit,*,iostat=ierr)LINE
        if (LINE.ne."z") then
          write(0,"(a,a)") "Error! Expected z, but read ", trim(LINE)
        end if
        backspace(bsunit)
        read(bsunit,*,iostat=ierr)LINE,k,rl, im
        if (k.ne.m) then
          write(0,"(a,i2,a,i2)") "Error. Expected z dimension ", m, "but got ", k
        end if
        bs(j)%z(m)=cmplx(rl,im,kind=8)
      end do
    end do
    
    close(bsunit)

    if (t==0.0d0) then
      if (method.eq."MCEv1") then
        do j=1,nbf
          do r=1,npes
            bs(j)%d_pes(r) = bs(j)%d_pes(r) * bs(j)%D_big
            bs(j)%a_pes(r) = bs(j)%d_pes(r) * cdexp(i*bs(j)%s_pes(r))
          end do
          bs(j)%D_big = (1.0d0,0.0d0)
        end do
      else
        do j=1,nbf
          dsum1 = (0.0d0,0.0d0)
          do r=1,npes
            dsum1 = dsum1 + bs(j)%d_pes(r)
            if (r.eq.in_pes) then
              bs(j)%d_pes(r) = (1.0d0,0.0d0)
            else
              bs(j)%d_pes(r) = (0.0d0,0.0d0)
            end if
            bs(j)%a_pes(r) = bs(j)%d_pes(r) * cdexp(i*bs(j)%s_pes(r))
          end do
          bs(j)%D_big = dsum1 * bs(j)%D_big
        end do       
      end if
    else
      do j=1,nbf
        if ((dble(bs(j)%D_big) /= 1.0d0).and.(method.eq."MCEv1")) then
          write(0,"(a)") "The D_big amplitudes are not compatible with MCEv1 propagation for t/=0"
          errorflag=1
          return
        else if ((dble(bs(j)%d_pes(1)) /= 1.0d0).and.(method.eq."CCS")) then
          write(0,"(a)") "The d_pes amplitudes are not compatible with CCS propagation for t/=0"
          errorflag=1
          return
        end if
      end do
    end if   
       
  end subroutine readbasis
  
!------------------------------------------------------------------------------------

  subroutine restartnum(rep,gen,restart)
  
    implicit none
    
    integer, intent (inout) :: restart
    integer, intent(in) :: rep
    character(LEN=1), intent(inout) :: gen

    real(kind=8) :: time    
    integer :: ierr,fileun
    character(LEN=100) :: LINE
    character(LEN=13)::filename
    logical :: file_exists
    
    if (errorflag.ne.0) return
    
    fileun=25433+rep
    write(filename, "(a,i3.3,a)") "Outbs-",rep,".out"
    inquire(file=filename,exist=file_exists)
    
    if (file_exists.eqv..false.) then
      gen="Y"
    else
      open(unit=fileun,file=filename,status="old",iostat=ierr)
      if (ierr/=0) then
        write (0,"(3a,i0)") "Error opening the ", trim(filename), " file. Ierr = ", ierr
        errorflag = 1
        return
      end if
      
      read(fileun,*,iostat=ierr)LINE

      do while ((LINE.ne."time").and.(ierr==0)) 
        read(fileun,*,iostat=ierr)LINE
      end do
      backspace(fileun)
      
      read(fileun,*,iostat=ierr) LINE, time
      
      if (time==timeend) then
        restart = 1
      else
        gen = "N"
      end if
      close(fileun)
    end if      
    
  end subroutine restartnum  
  
!--------------------------------------------------------------------------------------------------

  subroutine readtimepar   !   Level 1 Subroutine

    implicit none
    character(LEN=100)::LINE, LINE2
    integer::ierr, n

    if (errorflag .ne. 0) return

    ierr = 0
    n = 0

    OPEN(UNIT=135, FILE='prop.dat',STATUS='OLD', iostat=ierr)

    if (ierr .ne. 0) then
      write(0,"(a)") 'Error in opening prop.dat file'
      errorflag = 1
      return
    end if

    read(135,*,iostat=ierr)LINE

    do while (ierr==0)

      if(LINE== "dtmin") then
        backspace(135)
        read(135,*,iostat=ierr)LINE,dtmin 
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading minimum dt"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="dtmax") then
        backspace(135)
        read(135,*,iostat=ierr)LINE,dtmax
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading maximum dt"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="dtinit") then
        backspace(135)
        read(135,*,iostat=ierr)LINE,dtinit
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading initial dt"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="time_end") then
        backspace(135)
        read(135,*,iostat=ierr)LINE,timeend
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading end time for propagation"
          errorflag = 1
          return
        end if
        n = n+1    
      else if (LINE=="time_start") then
        backspace(135)
        read(135,*,iostat=ierr)LINE,timestrt
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading starting time for propagation"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="step") then
        backspace(135)
        read(135,*,iostat=ierr)LINE,LINE2
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading step type for propagation"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'S').or.(LINE2(1:1).eq.'s')) then
          step = "S"
        else if ((LINE2(1:1).eq.'A').or.(LINE2(1:1).eq.'a')) then
          step = "A"
        else
          write(0,"(a,a)") "Error reading step type. Expected 'Static' or 'Adaptive', but got ", trim(LINE2)
          errorflag = 1
        end if
        n = n+1
      end if    

      read(135,*,iostat=ierr)LINE

    end do

    if (n.ne.6) then
      write(0,"(a)") "Not all required variables read in readtimepar subroutine."
      errorflag = 1
      return
    end if

    return

  end subroutine readtimepar
  
!--------------------------------------------------------------------------------------------------

  subroutine readclone(clonenum, reps, clone)
  
    implicit none
    
    integer, dimension(:), intent(inout) :: clonenum, clone
    integer, intent(in) :: reps
    
    integer :: j, k, p, q, n, ierr
    character(LEN=255) :: filename
    character(LEN=3) :: rep
    
    if (errorflag .ne. 0) return
    
    write(rep,"(i3.3)") reps
    filename="clonearr-"//rep//".out"
    
    n = 90254+reps
    
    open(unit=n,file=trim(filename),status="old",iostat=ierr)
    if (ierr/=0) then
      write (0,"(3a,i0)") "Error opening ", trim(filename), " file. Ierr was ", ierr
      errorflag = 1
      return
    end if
    
    do j=1,size(clonenum)
      read (n,*,iostat=ierr) k, p, q
      if (ierr/=0) then
        write (0,"(3a,i0,a,i0)") "Error reading from ", trim(filename), " file when j = ", j, &
                                  ". Ierr was ", ierr
        errorflag = 1
        return
      end if
      if (k.ne.j) then
        write (0,"(3a,i0,a,i0)") "Error reading from ", trim(filename), " file. Expected j = ", j, &
                                  " but read k = ", k
        errorflag = 1
        return
      end if  
      clonenum(j) = p
      clone(j) = q
    end do
    
    close(n)
        
    return
  
  end subroutine readclone  
    
!*************************************************************************************************!

END MODULE readpars
