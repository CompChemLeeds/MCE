module test1

contains

function arraytst(a,b,c)

  implicit none
  
  real(kind=8), dimension(:,:,:), allocatable :: arraytst
  integer, intent(in) :: a,b,c
  
  allocate(arraytst(a,b,c))
  
  arraytst = (dble(a)/dble(b))*dble(c)
  
  return
  
end function arraytst

end module test1

Program Testsuite

  use test1

  implicit none
  
  real(kind=8) :: time, testsum
  real(kind=8), dimension(:,:,:), allocatable :: test
  character(LEN=10) :: timestr
  integer :: i, j, k
    real :: r
    integer(kind=8)::ranseed
    integer,allocatable::ranseed_array(:)
    integer:: ranseed_size, rint
    integer::clock, m, ierr
    character(LEN=100)::LINE, command

    call random_seed(size=ranseed_size)
    allocate(ranseed_array(ranseed_size))
    call system_clock(count=clock)
    ranseed_array = clock + 37* (/ (i-1,i=1,ranseed_size) /)
    call random_seed(put=ranseed_array)
    deallocate(ranseed_array)
    CALL RANDOM_NUMBER(r)
    r=r*3232768.0
    rint = int(r)

  time = 4.99750040d0
  
  i=2
  j=3
  k=4
  
    write(LINE,"(i3.3)") j*k 
    call system("cp ./sb.f90 ./sb"//trim(LINE)//".f90")

  allocate(test(i,j,k))

  test=arraytst(i,j,k)
  
  testsum = sum(test(i,j,:))
  
  write(6,*) test
  
  write(6,*) testsum

  write(timestr,"(i3.3,f0.4)") int(time), time-int(time)
  write(6,"(3a)") "___",timestr,"___"
  write(6,"(3a)") "___",trim(timestr),"___"
 
  stop

end program Testsuite
