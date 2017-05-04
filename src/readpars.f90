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
!*      1) Reading the run conditions (debug,gen,prop,cmprss,method,reptot,conjflg)
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

      if(LINE=='debug') then
        backspace(140)
        read(140,*,iostat=ierr)LINE,debug
        if (ierr.ne.0) then
          write(0,"(a)") "Error reading debug status"
          errorflag = 1
          return
        end if
        if ((debug.ne.0).and.(debug.ne.1)) then
          write(0,"(a)") "Error in debug state. Debug value should be only 0 (for off) or 1 (for on)"
          write(0,"(a,i3)") "Debug value read was ", debug
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=='gen') then
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
        else if ((LINE2(1:9).eq.'aimc-mce1').or.(LINE2(1:9).eq.'AIMC-MCE1')) then
          method = "AIMC1"
        else if ((LINE2(1:9).eq.'aimc-mce2').or.(LINE2(1:9).eq.'AIMC-MCE2')) then
          method = "AIMC2"
        else
          write(0,"(a,a,a)") "Error. Method must be MCEv1, MCEv2, AIMC-MCE1, ",&
                        "AIMC-MCE2 or CCS. Read ", trim(LINE2)
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

    if (n.ne.8) then
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
      end if
      read(127,*,iostat=ierr) LINE

    end do

    close(127)

    if (n.ne.1) then
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
      else if (LINE=="gridsp") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,initsp
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid spacing"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="psizex") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,psizex
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid size in p coordinate"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="qsizex") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,qsizex
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid size in q coordinate"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="psizey") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,psizey
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid size in p coordinate"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="qsizey") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,qsizey
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid size in q coordinate"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="psizez") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,psizez
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid size in p coordinate"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="qsizez") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,qsizez
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading grid size in q coordinate"
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
      else if (LINE=="def_stp") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,def_stp
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading default number of basis functions per train"
          errorflag = 1
          return
        end if
        n = n+1
      else if (LINE=="nbfadapt") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,LINE2
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading bf adapt flag"
          errorflag = 1
          return
        end if
        if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
          nbfadapt = "YES"
        else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
          nbfadapt = "NO"
        else
          write(0,"(a,a)") "Error. nbfadapt value must be YES/NO. Read ", trim(LINE2)
        end if
        n = n+1
      else if (LINE=="nbfepsilon") then
        backspace(127)
        read(127,*,iostat=ierr)LINE,bfeps
        if(ierr.ne.0) then
          write(0,"(a)")  "Error reading bf adapt cutoff parameter"
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
          cloneflg = "BLIND"
        else
          write(0,"(a,a)") "Error. cloneflg value must be YES/NO/BLIND. Read ", trim(LINE2)
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
    
    if ((basis.ne.'TRAIN').and.(basis.ne.'train').and.(basis.ne.'SWARM').and.(basis.ne.'swarm').and.(basis.ne.'SWTRN')&
      .and.(basis.ne.'swtrn').and.(basis.ne.'GRID').and.(basis.ne.'grid').and.(basis.ne.'GRSWM').and.(basis.ne.'grswm')) then
      write(0,"(a,a)") "Invalid value for basis. Must be TRAIN/SWARM/GRID/SWTRN/GRSWM and all upper/lower case Value is ", basis
      errorflag = 1
      return
    else if ((basis.eq.'TRAIN').or.(basis.eq.'train'))then
      basis = 'TRAIN'
    else if ((basis.eq.'SWARM').or.(basis.eq.'swarm'))then
      basis = 'SWARM'
    else if ((basis.eq.'GRID').or.(basis.eq.'grid'))then
      basis = 'GRID'
    else if ((basis.eq.'SWTRN').or.(basis.eq.'swtrn'))then
      basis = 'SWTRN'
    else if ((basis.eq.'GRSWM').or.(basis.eq.'grswm'))then
      basis = 'GRSWM'
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
    
    if (n.ne.20) then
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
    
    if (basis.eq."GRID") then
      if (initsp .lt. 0.8d0) then
        write(0,"(a)") "Error! Grid points are too close together"
        errorflag=1
        return
      else if (initsp .gt. 1.85d0) then
        write(0,"(a)") "Error! Grid points are too widely spaced"
      end if
      if ((ndim.ne.1).and.(ndim.ne.3)) then
        write(0,"(a)") "ndim is neither 1 nor 3. This is currently invalid."
        errorflag=1
        return
      else if (ndim==1) then
        if ((qsizex.gt.qsizey).and.(qsizex.gt.qsizez).and.(psizex.gt.psizey).and.(psizex.gt.psizez)) then
          qsizey=0
          qsizez=0
          psizey=0
          psizez=0
        else if ((qsizey.gt.qsizex).and.(qsizey.gt.qsizez).and.(psizey.gt.psizex).and.(psizey.gt.psizez)) then
          qsizex=0
          qsizez=0
          psizex=0
          psizez=0
        else if ((qsizez.gt.qsizex).and.(qsizez.gt.qsizey).and.(psizez.gt.psizex).and.(psizez.gt.psizey)) then
          qsizex=0
          qsizey=0
          psizex=0
          psizey=0
        else 
          write(0,"(a)") "Error! The largest grid dimensions do not match. Check the values."
          errorflag = 1
          return
        end if
      end if
      if ((mod(qsizez,2).ne.mod(psizez,2)).and.(mod(qsizey,2).ne.mod(psizey,2)).and.(mod(qsizex,2).ne.mod(psizex,2))) then
        write(0,"(a)") "Error! Grid is not symmetrically spaced around initial CS."
        write(0,"(a)") "There should be equal distance between initial CS and the four closest grid points"
        write(0,"(a)") "This equates to  all psize and qsize values being either all even or all odd."
        errorflag=1
        return
      end if
      if ((max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)).ne.in_nbf) then
        write(6,"(a)") "Warning! Grid size does not match in_nbf. in_nbf should be product of all psize and qsize"
        if (mod(qsizez,2)==1) then
          in_nbf = (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1))
          write(6,"(a,i0,a)") "Altering in_nbf to ", &
               (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)), "..."
        else
          in_nbf = (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)) + 1
          write(6,"(a,i0,a)") "Altering in_nbf to ", &
                  (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)) + 1, "..."
        end if
      end if
    end if
    
    if (basis.eq."GRSWM") then
      if (initsp .lt. 0.8d0) then
        write(0,"(a)") "Error! Grid points are too close together"
        errorflag=1
        return
      end if
      if (ndim.ne.3) then
        write(0,"(a)") "ndim is not 3. This is currently invalid."
        errorflag=1
        return
      end if
      if (mod(qsizez,2).ne.mod(psizez,2)) then
        write(0,"(a)") "Error! Grid is not symmetrically spaced around initial CS."
        write(0,"(a,a)") "There should be equal distance between initial CS and the four",&
                           " closest grid points in the grid dimension (z)"
        write(0,"(a)") "This equates to both psizez and qsizez values being even or odd."
        errorflag=1
        return
      end if
      if ((qsizez*psizez).ne.in_nbf) then
        write(6,"(a)") "Warning! Grid size does not match in_nbf. in_nbf should be product of psizez and qsizez"
        if (mod(qsizez,2)==1) then
          in_nbf = qsizez*psizez
          write(6,"(a,i0,a)") "Altering in_nbf to ", qsizez*psizez, "..."
        else
          in_nbf = qsizez*psizez + 1
          write(6,"(a,i0,a)") "Altering in_nbf to ", qsizez*psizez+1, "..."
        end if
      end if
    end if

    if ((basis.eq."TRAIN").or.(basis.eq."SWTRN")) then
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
      if ((method.ne."MCEv2").and.(method.ne."CCS").and.(method.ne."AIMC2")) then
        write(0,"(a)") "Error! Trains can only work with MCEv2, AIMC-MCE (second pass) or CCS."
        write(0,"(a)") "MCEv1 and AIMC-MCE (first pass) cannot calculate the amplitudes correctly"
        errorflag = 1
        return
      end if        
    end if
    
    !!!!! Method dependent checks !!!!!!!
    
    if (method.eq."AIMC1") then
      if (basis.ne."SWARM") then
        write(0,"(a)") "Error! The AIMC-MCE first pass should only be carried out using a swarm basis set"
        errorflag = 1
        return
      end if
      if (prop == "N") then
        write (0,"(a)") "Error! The AIMC-MCE first pass should be propagated. Propagation flag set to N"
        errorflag = 1
        return
      end if
      if (cloneflg == "NO") then
        write (0,"(a)") "AIMC-MCE first pass requires cloning to be enabled."
        write (0,"(a)") "Enabling now"
        cloneflg = "YES"
      end if
      if (mod(def_stp,2)==0) then
        write(0,"(2(a,i0))") "For AIMC-MCE, def_stp should be odd. Changed from ", def_stp, " to ", def_stp+1
        def_stp = def_stp+1
      end if
    end if
    
    if (method.eq."AIMC2") then
      if (basis.ne."SWTRN") then
        write (0,"(a)") "Error! Basis must be a swarm of trains for AIMC-MCE second pass"
        errorflag = 1
        return
      end if
      if (prop.eq."N") then
        write (0,"(a)") "Error! The AIMC-MCE second pass should be propagated. Propagation flag set to N"
        errorflag = 1
        return
      end if
      if (mod(def_stp,2)==0) then
        write(0,"(2(a,i0))") "For AIMC-MCE, def_stp should be odd. Changed from ", def_stp, " to ", def_stp+1
        def_stp = def_stp+1
      end if
!      if (gen.eq."Y") then
!        write (0,"(a,a)") "Error! The AIMC-MCE second pass relies on precalculated basis functions, but",&
!                        " generation flag set to Y"
!        errorflag = 1
!        return
!      end if
    end if
    
    !!!!!! Basis set change parameter check!!!!!!

    if (cloneflg.ne.'NO') then
      if ((method.ne."MCEv2").and.(method.ne."AIMC1")) then
        write (0,"(a)") "Cloning can only work with MCEv2 or AIMC-MCE (first pass)"
        errorflag = 1
        return
      else if (method.eq."MCEv2") then
        if ((basis.ne."SWARM").and.(basis.ne."SWTRN").and.(basis.ne."TRAIN")) then
          write (0,"(a)") "Cloning with MCEv2 can only work on a swarm, a train or a swarm of trains"
          errorflag = 1
          return
        end if
      else
        if (basis.ne."SWARM") then
          write(0,"(a)") "AIMC-MCE (first pass) should be done only with a swarm of train centres."
          errorflag = 1
          return
        end if
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
        
    if (nbfadapt.eq."YES") then
      if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
        write(0,"(a)") "Adaptive basis set chosen but gridding disabled. This is currently an invalid combination."
        write(0,"(a,a)") "basis should be GRID or GRSWM. Value was ", basis 
        errorflag=1
        return
      end if
      if (bfeps.lt.0.0d0) then
        write(0,"(a)") "Basis set adaptive cutoff parameter less than zero. This is not valid"
        errorflag=1
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
        if ((method.ne."MCEv1").and.(method.ne."MCEv2").and.(method.ne."AIMC1").and.(method.ne."AIMC2")) then
          write(0,"(a)") "Spin Boson model can only be simulated by MCEv1 or MCEv2, or through AIMC-MCE"
          errorflag = 1
          return
        end if
        if (basis.eq."GRID") then
          write(0,"(a)") "This method must not use a static grid."
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
      case ("FP")
        if (npes.ne.1) then
          write(0,"(a)") "Free Particle only valid for 1 PES"
          errorflag = 1
          return
        end if
        if (trim(method).ne."CCS") then
          write(0,"(a)") "Free Particle can only be simulated by CCS"
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
      case ("IV")
        if (npes.ne.1) then
          write(0,"(a)") "Inverted Gaussian only valid for 1 PES"
          errorflag = 1
          return
        end if
        if ((ndim.ne.1).and.(ndim.ne.3)) then
          write(0,"(a)") "Inverted Gaussian is only valid for 1 or 3 dimensional"
          errorflag = 1
          return 
        end if
        if (trim(method).ne."CCS") then
          write(0,"(a)") "Inverted Gaussian can only be simulated by CCS"
          errorflag = 1
          return  
        end if 
        if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
          write(0,"(a)") "This method must use a static grid."
          errorflag = 1
          return
        end if
      case ("CP")
        if (npes.ne.1) then
          write(0,"(a)") "Coulomb Potential only valid for 1 PES"
          errorflag = 1
          return
        end if
        if (mod(qsizez,2)==1) then
          write(0,"(a)") "Odd parity for the grid parameters will result in a point on the singularity"
          write(0,"(a)") "This would make reprojection impossible, as well as seriously affecting results"
          write(0,"(a)") "Change the grid parameters and restart."
          errorflag=1
          return
        end if
        if (in_nbf.eq.(qsizez*psizez+1)) then
          write(6,"(a)") "The in_nbf value has been reset so that there is a point on the singularity."
          write(6,"(a)") "This will make reprojection impossible. Resetting to even parity."
          in_nbf = qsizez*psizez
        end if
        if (ndim.ne.3) then
          write(0,"(a)") "Coulomb Potential is only valid in 3 dimensions"
          errorflag = 1
          return 
        end if
        if (trim(method).ne."CCS") then
          write(0,"(a)") "Coulomb Potential can only be simulated by CCS"
          errorflag = 1
          return  
        end if 
        if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
          write(0,"(a)") "This method must use a static grid."
          errorflag = 1
          return
        end if
      case ("HH")
        if (npes.ne.1) then
          write(0,"(a)") "Inverted Gaussian only valid for 1 PES"
          errorflag = 1
          return
        end if
        if (trim(method).ne."CCS") then
          write(0,"(a)") "Henon-Heiles Potential can only be simulated by CCS"
          errorflag = 1
          return  
        end if 
        if ((ndim.ne.2).and.(ndim.ne.6).and.(ndim.ne.10)) then
          write(0,"(a)") "Henon-Heiles potential is only valid currently for the 2,6,or 10 dimensional systems"
          errorflag = 1
          return 
        end if
        if (basis.eq."GRID") then
          write(0,"(a)") "This method must not use a static grid."
          errorflag = 1
          return
        end if 
      case default
        write(0,"(a)") "System is not recognised. Value is ", sys
        errorflag = 1
        return
    end select 
    
    initsp=initsp*dsqrt(2.0d0) 
    
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
      end if

      read(127,*,iostat=ierr) LINE

    end do

    close(127)

    sigp = sqrt(1.0d0/(2.0d0*gam))
    sigq = sqrt(gam/2.0d0)

    if (n.ne.3) then
      write(0,"(a)") "Not all required variables read in readzparams subroutine."
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
      write(0,"(a,i2,a)") "Error in reading parameters. Only ", n, " of 6 parameters read."
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

  subroutine constrtrain(bs, x, time, reps, mup, muq, rkstp, genflg,nbf,map_bfs)
  
    implicit none
    type(basisfn), intent(inout), dimension (:), allocatable :: bs
    real(kind=8),  intent(inout), dimension (:), allocatable :: mup, muq
    real(kind=8),  intent(inout) :: time
    integer, intent(inout), dimension (:,:) :: map_bfs
    integer, intent(in) :: x, reps,rkstp, genflg, nbf
    
    real(kind=8) :: t, rl, im
    integer, dimension (:), allocatable :: bfs
    integer :: carriage, stpback, bsunit, j, k, l, m, n, r, p, q, ierr, nbftrk, nbftrk2, temppar, maxbf
    character(LEN=255) :: LINE
    character(LEN=21) filename
    character(LEN=3) :: rep
    character(LEN=5) :: step, tempchar
    character(LEN=1) :: rkstp_char
    
    if (errorflag.ne.0) return

    stpback = ((def_stp-1)/2)*trainsp
    
    carriage = x+stpback
    write(rep,"(i3.3)") reps
    write(rkstp_char,"(i1.1)") rkstp
    
    if (allocated(bs)) then
      write (0,"(a)") "The basis set in the train construction subroutine should be unallocated at call time"
      write (0,"(a)") "This means either a dummy basis set variable or the unallocated initial basis should "
      write (0,"(a)") "be sent to the subroutine, as it will be allocated based on the input file sizes"
      errorflag = 1
      return
    end if
    
    do j=1,def_stp
      bsunit=(7065+carriage)*reps
      if (carriage.lt.0) then
        write(step,"(i5.4)") carriage
      else
        write(step,"(i5.5)") carriage
      end if
      filename="Outbs-"//trim(rep)//"-"//trim(step)//"-"//trim(rkstp_char)//".out"
      open(unit=bsunit,file=trim(filename),status="old",iostat=ierr)
      if (ierr/=0) then
        write(0,"(a,a,a,i0)") "Error opening file ", trim(filename), " in step ", x
        write(0,"(a,i0)") "ierr was ", ierr 
        errorflag = 1
        return
      end if
      carriage = carriage - trainsp
    end do

    if ((x==1).and.(rkstp==0).and.(genflg==1)) then 
    
      bsunit = 7065*reps+1
      n=0
      
      read(bsunit,*,iostat=ierr)LINE

      do while ((LINE.ne."zinit").and.(ierr==0))
        if (LINE=="ndof") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,ndim
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading ndim"
            errorflag = 1
            return
          end if
          write(6,"(a,i0)") "ndim   = ", ndim
          n = n+1
        else if (LINE=="nconf") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,npes
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading npes"
            errorflag = 1
            return
          end if
          write(6,"(a,i0)") "npes   = ", npes
          n = n+1
        else if (LINE=="initial_PES") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,in_pes
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading in_PES"
            errorflag = 1
            return
          end if
          write(6,"(a,i0)") "in_pes = ", in_pes
          n = n+1
        else if (LINE=="nbasisfns") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,in_nbf
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading in_nbf"
            errorflag = 1
            return
          end if
          write(6,"(a,i0)") "in_nbf = ", in_nbf
          n = n+1
        else if (LINE=="matfun") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,matfun
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading matrix function"
            errorflag = 1
            return
          end if
          write(6,"(a,a)") "matfun = ", matfun
          n = n+1
        else if (LINE=="time") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,time
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading time"
            errorflag = 1
            return
          end if
          write(6,"(a,es16.8e3)") "time   = ", time
          n = n+1
        end if
        read(bsunit,*,iostat=ierr)LINE
      end do 

      if (n.ne.6) then
        write(0,"(a,i0,a)") "Error in reading parameters from initial point. ", n, " of 6 parameters read."
        errorflag = 1
        return
      end if 
      
      do j=1,def_stp
        do k=1,in_nbf
          map_bfs(j,k) = ((j-1)*in_nbf)+k
        end do
      end do
      
      allocate (mup(ndim), stat=ierr)
      if (ierr == 0) allocate (muq(ndim), stat=ierr)
      if (ierr/=0) then
        write(0,"(a)") "Error in allocation of mup and muq"
        errorflag=1
      end if
      
      backspace(bsunit)
      
      do m=1,ndim
        read(bsunit,*,iostat=ierr)LINE, j, muq(j), mup(j)
        if(ierr.ne.0) then
          write(0,"(a,i0)")  "Error reading zinit value ", m
          errorflag = 1
          return
        end if
        if(m.ne.j) then
          write(0,"(2(a,i0))")  "Error! Count mismatch in zinit. Expected ", m, " but read ", j
          errorflag = 1
          return
        end if
      end do
      
      rewind (bsunit)    
      
    end if          !!!!!! End of initial reading of basis set parameters
    
    carriage = x+stpback
    nbftrk = 0
    allocate(bfs(def_stp), stat = ierr)
    if (ierr/=0) then
      write (0,"(a,i0)") "Error allocating the basis functions. ierr was ", ierr
      errorflag = 1
      return
    end if
    
    do j=1,def_stp

      bsunit=(7065+carriage)*reps
      n=0
     
      read(bsunit,*,iostat=ierr)LINE

      do while ((LINE.ne."basis").and.(ierr==0))
        if (LINE=="ndof") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,temppar
          if(ierr.ne.0) then
            write(0,"(a,i0)")  "Error reading ndim for carriage ", carriage
            errorflag = 1
            return
          end if
          if (temppar.ne.ndim) then
            write(0,"(2(a,i0))") "Mismatch in ndim values! Expected ", ndim, " but got ", temppar
            errorflag = 1
            return
          end if
          n = n+1
        else if (LINE=="nconf") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,temppar
          if(ierr.ne.0) then
            write(0,"(a,i0)")  "Error reading npes for carriage ", carriage
            errorflag = 1
            return
          end if
          if (temppar.ne.npes) then
            write(0,"(2(a,i0))") "Mismatch in npes values! Expected ", npes, " but got ", temppar
            errorflag = 1
            return
          end if
          n = n+1
        else if (LINE=="nbasisfns") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,temppar
          if(ierr.ne.0) then
            write(0,"(a,i0)")  "Error reading swarm nbf for carriage ", carriage
            errorflag = 1
            return
          end if
          nbftrk = nbftrk + temppar
          bfs(j) = temppar
          n = n+1
        else if (LINE=="initial_PES") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,temppar
          if(ierr.ne.0) then
            write(0,"(a)")  "Error reading in_pes for carriage ", carriage
            errorflag = 1
            return
          end if
          if (temppar.ne.in_pes) then
            write(0,"(a)") "Mismatch in in_pes values! Expected ", in_pes, " but got ", temppar
            errorflag = 1
            return
          end if
          n = n+1
        else if (LINE=="matfun") then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,tempchar
          if(ierr.ne.0) then
            write(0,"(a,i0)")  "Error reading matrix function for carriage ", carriage
            errorflag = 1
            return
          end if
          if (tempchar.ne.matfun) then
            write(0,"(4a)") "Mismatch in matfun values! Expected ", matfun, " but got ", tempchar
            errorflag = 1
            return
          end if
          n = n+1
        else if ((LINE=="time").and.(carriage.eq.x)) then
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,t
          if(ierr.ne.0) then
            write(0,"(a,i0)")  "Error reading time in carriage ", carriage
            errorflag = 1
            return
          end if
          if (t-time.gt.dtinit/1.0d3) then
            write(0,"(2(a,es16.8e3))") "Mismatch in time values! Expected ", time, " but got ", t
            write(0,"(3(a,i0))") "This happened in rep ", reps, " at step ", x, " and rkstp ", rkstp
            errorflag = 1
            return
          end if
          n = n+1
        end if
        read(bsunit,*,iostat=ierr)LINE
      end do

      if ((n.ne.5).and.(carriage.ne.x)) then
        write(0,"(a,i0,a)") "Error in reading parameters. Only ", n, " of 5 parameters read."
        errorflag = 1
        return
      else if ((n.ne.6).and.(carriage.eq.x)) then
        write(0,"(a,i0,a)") "Error in reading parameters. Only ", n, " of 6 parameters read."
        errorflag = 1
        return
      end if 
      backspace (bsunit)
      carriage = carriage - trainsp
    end do      

    call allocbs(bs,nbftrk)
    
    nbftrk2 = 0
    carriage = x+stpback
    p=0
    n=0
    maxbf = maxval(map_bfs)
 
    do l=1,def_stp
      bsunit=(7065+carriage)*reps
      do j=1,bfs(l)
        n=n+1
        if (map_bfs(l,j)==0) then
          p=p+1
          q=maxbf+p
          map_bfs(l,j) = q
        else
          q=map_bfs(l,j)
        end if  
        read(bsunit,*,iostat=ierr)LINE,k
        if (ierr/=0) then
          write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
          write(0,"(a)") "Was trying to read basis number line" 
          write (0,"(a,i0)") "ierr was ", ierr
          errorflag = 1
          return
        end if
        if(k.ne.j) then
          write(0,"(a,i0,a,i0)") "Error. Expected basis function ", j, " but got ", k
          errorflag = 1
          return
        end if
        if (LINE.ne."basis") then
          write(0,"(a,a)") 'Error. Expected to read "basis ', j, '" but read ',trim(LINE), ' ', k 
          errorflag = 1
          return
        end if
        read (bsunit,*,iostat=ierr)LINE
        if (ierr/=0) then
          write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
          write(0,"(a)") "Was trying to D line" 
          write (0,"(a,i0)") "ierr was ", ierr
          errorflag = 1
          return
        end if
        if (LINE.ne."D") then
          write(0,"(a,a)") "Error! Expected D but read ", trim(LINE)
        end if
        do r=1,npes
          read(bsunit,*,iostat=ierr)LINE
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read a_pes identifier for pes ", r 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (LINE.ne."a") then
            write(0,"(a,a)") "Error! Expected a but read ", trim(LINE)
          end if
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,k,rl, im
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read a_pes line for pes ", r 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (k.ne.r) then
            write (0,"(a,i2,a,i2)") "Error. Expected a from pes ", r, "but got ", k
          end if
          bs(q)%a_pes(r)=cmplx(rl,im,kind=8)
        end do
        do r=1,npes
          read(bsunit,*,iostat=ierr)LINE
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read d_pes identifier for pes ", r 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (LINE.ne."d") then
            write(0,"(a,a)") "Error! Expected d but read ", trim(LINE)
          end if
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,k,rl, im
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read d_pes line for pes ", r 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (k.ne.r) then
            write(0,"(a,i2,a,i2)") "Error. Expected d from pes ", r, "but got ", k
          end if
          bs(q)%d_pes(r)=cmplx(rl,im,kind=8)
        end do
        do r=1,npes
          read(bsunit,*,iostat=ierr)LINE
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read s_pes identifier for pes ", r 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (LINE.ne."s") then
            write(0,"(a,a)") "Error! Expected s but read ", trim(LINE)
          end if
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,k,rl
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read s_pes line for pes ", r 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (k.ne.r) then
            write(0,"(a,i2,a,i2)") "Error. Expected s from pes ", r, "but got ", k
          end if
          bs(q)%s_pes(r)=rl
        end do
        do m=1,ndim
          read(bsunit,*,iostat=ierr)LINE
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read z identifier for dof ", m 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (LINE.ne."z") then
            write(0,"(a,a)") "Error! Expected z, but read ", trim(LINE)
          end if
          backspace(bsunit)
          read(bsunit,*,iostat=ierr)LINE,k,rl, im
          if (ierr/=0) then
            write(0,"(a,i0,a,i0)") "Error reading Outbs file ", carriage, " at basis number ", j
            write(0,"(a,i0)") "Was trying to read z line for dof ", m 
            write (0,"(a,i0)") "ierr was ", ierr
            errorflag = 1
            return
          end if
          if (k.ne.m) then
            write(0,"(a,i2,a,i2)") "Error. Expected z dimension ", m, "but got ", k
          end if
          bs(q)%z(m)=cmplx(rl,im,kind=8)
        end do
      end do

      carriage = carriage-trainsp
      close(bsunit)
    end do 
    
    if (maxval(map_bfs).ne.maxbf+p) then
      write(0,"(a)") "The two nbf trackers do not agree. Problem with cloning indices!"
      write(0,"(2(a,i0))") "maxval(map_bfs) was ", maxval(map_bfs), " and maxbf+p was ", maxbf+p 
      errorflag = 1
      return
    end if  

    return
    
  end subroutine constrtrain

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
