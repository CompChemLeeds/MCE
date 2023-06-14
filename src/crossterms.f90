!!****************************************************************
! First attempt at implementing crossterms, ran into issues with makefiles and passing back and forth between python and 
! fortran. Archieved as maybe some useful code written here to build off of later.
! 14/10/2022 - RB
!! *******************************************************************


Program crossterms
    use bsetgen    ! basis set generation module
    use bsetalter  ! module to change the size/position of the basis set
    use Ham        ! General hamiltonian module, inc. overlap calculations etc
    use globvars   ! Global variables, used in all modules
    use readpars   ! Module to read parameters from input file
    use outputs    ! Module to output data
    use alarrays   ! Array allocation module for defined types
    use Chks       ! A set of checks which ensure the program is running correctly
    use propMCE    ! The time propagation controller
    use redirect   ! Module which directs functions to the system specific variants

    implicit none

    type(basisfn), dimension (:), allocatable :: bs1, bs2
    real(kind=8) ::normw1, normw2
    real(kind=8) ::  crossterm1, crossterm2
    integer :: j, k, ierr, nbf
    complex(kind=8) :: sumamps1, sumamps2, ovrlp
    character(len=100) :: filename1, filename2
    character(len=3) :: nbfin
    character(len=10) :: normwin1, normwin2
    

    ! need to write the inputs for the names of the two files as well as the normweights

    call getarg(0, filename1)
    if (ierr==-1) then
      write (0,"(a,a)") "Error! Could not read first argument for ", trim(filename1)
      write (0,"(a)") "This should be the name of a basis set tracker."
      stop
    end if
    

    call getarg(1, normwin1)
    if (ierr==-1) then
      write (0,"(a,a)") "Error! Could not read first argument for ", normw1
      write (0,"(a)") "This should be the first normweighting."
      stop
    end if
    write(normwin1, *) normw1


    call getarg(2, filename2)
    if (ierr==-1) then
      write (0,"(a,a)") "Error! Could not read first argument for ", trim(filename2)
      write (0,"(a)") "This should be the name of a basis set tracker."
      stop
    end if

    call getarg(3, normwin2)
    if (ierr==-1) then
      write (0,"(a,a)") "Error! Could not read first argument for ", normw2
      write (0,"(a)") "This should be the second normweighting."
      stop
    end if
    write(normwin2, *) normw2
    call getarg(4,nbfin)
    if (ierr==-1) then
      write (0,"(a,a)") "Error! Could not read first argument for ", nbfin
      write (0,"(a)") "This should be the number of basis functions."
      stop
    end if
    write(nbfin, *) nbf



    sumamps1 = (0.0d0, 0.0d0)
    sumamps2 = (0.0d0, 0.0d0)

    ! Now to read in the basis sets themselves.
    call allocbs(bs1,nbf)
    call readcontbasis(bs1,trim(filename1),nbf)

    call allocbs(bs2,nbf)
    call readcontbasis(bs2,trim(filename2), nbf)

    !Calculating the overlap between the basis sets and summing over amplitudes

    do j=1,nbf
      do k=1,nbf
        ovrlp = ovrlp + ovrlpij(bs1(k)%z(:), bs2(j)%z(:))
        sumamps1 = sumamps1 + (dconjg(bs1(j)%a_pes(1))*bs2(k)%a_pes(1))
        sumamps2 = sumamps2 + (dconjg(bs1(j)%a_pes(2))*bs2(k)%a_pes(2))
      end do
    end do

    ! Calculating the cross terms

    crossterm1 = 2*real(sqrt((1-normw1)*normw2)*sumamps1*ovrlp)
    crossterm2 = 2*real(sqrt(normw1*(1-normw2))*sumamps2*ovrlp)
    
end program 
    
