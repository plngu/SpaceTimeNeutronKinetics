program SpaceTimeKinetic

  real           :: start, finish, a(2), seconds
  integer        :: inp, channel, channels(25), numErr
  integer        :: hours, minutes
  character(10)  :: date, time, zone

  call cpu_time(start)
  inp = 0

  call b08ini(inp)

  if (inp.lt.0) then
    print 1, inp
    stop
  endif

  channels(:)  = 0
  channels(1)  = 1
  channels(17) = 17
  channels(21) = 21

  do channel = 1, 25
    if (channels(channel) == 0) cycle
    call b08(0,channels(channel),numErr)
    if (numErr < 0) then
      print 2, numErr, channel
      stop
    endif
  enddo

  call main
  call rdind2(17,'GGRR',a,2,1,0)

  do channel = 1, 25
    if (channels(channel) == 0) cycle
    call b08cls(0,channels(channel),numErr)
    if (numErr < 0) then
      print 4, numErr
      stop
    endif
  enddo

  call b08hlt(numErr)

  call cpu_time(finish)
  call date_and_time(date,time,zone)

  print 5, date(7:8), date(5:6),date(1:4),time(1:2),time(3:4),time(5:)

  hours   = int((finish-start)/3600)
  minutes = int((finish-start)/60)-hours*60
  seconds = (finish-start)-minutes*60-hours*3600

  print 6, hours, minutes, seconds

1 format(' B20 Error in: "b08ini", numErr = ',i8)
2 format(' B20 Error in: "b08", numErr = ',i8,' channel=',i8)
3 format(/' Calculation time "STK" = ',f8.3,' seconds')
4 format(' B20 Error in: "b08cls", numErr = ',i8)
5 format(' Day ',a,' month ',a,' year ',a,' | time ',a,' hours ',     &
         a,' minutes ',a,' seconds')
6 format(' Calculation time "STK" = ',i4,' hours ',                   &
         i2,' minutes ',f8.3,' seconds')

end
