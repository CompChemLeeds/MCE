Program Testsuite

  implicit none
  
  real(kind=8) :: time
  character(LEN=10) :: timestr

  time = 4.99750040d0

  write(timestr,"(i3.3,f0.4)") int(time), time-int(time)
  write(6,"(3a)") "___",timestr,"___"
  write(6,"(3a)") "___",trim(timestr),"___"
 
  stop

end program Testsuite
