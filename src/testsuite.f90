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

  time = 4.99750040d0
  
  i=2
  j=3
  k=4

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
