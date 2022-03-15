subroutine main
  call input
  call getDataFromXT26
  call timeFinder
  call initSource
  call findSolve
  return
end



subroutine getDataFromXT26

  use STK

  numMatr    = params(1)                                                ! кол-во матриц на интревал
  minTime    = time(1)                                                  ! минимальное время жизни в задаче
  maxTime    = time(2)                                                  ! максимальное время жизни в задаче
  linTime    = 0.0d0                                                    ! начало расчёта
  dt         = time(3)                                                  ! линейный шаг по времени равен минимальному времени жизни
  stopTime   = time(4)
  steps      = nint(maxTime/dt)                                         ! число линейных шагов по времени
  numOfZones = params(2)                                                ! количество зон в модели xt26
  modStep    = params(3)
  neutronGenerationCount = params(4)

  allocate (lnTime(numMatr),                                          &
            initAbsorbMatrix(numOfZones,numOfZones,0:numMatr),        &
            initFissionMatrix(numOfZones,numOfZones,0:numMatr),       &
            leakageMatrix(numOfZones,numMatr),                        &
            absMatrix(numOfZones), fisMatrix(numOfZones),             &
            leakMatrix(numOfZones),                                   &
            typeSource(numOfZones),                                   &
            typeSourceD(numOfZones))

  print 7
!$ Начальный источник в соответствии с плотностью деления
  typeSource  = 0.0
  typeSourceD = 0.0d0
  call rdind2(17,'GGRR',typeSource,numOfZones,1,0)
  typeSourceD = typeSource(1:numOfZones)

  print*, 'Init fission source'
  print 6, typeSourceD(1:numOfZones)
  noneZero = countNonZeroVectorElem(typeSourceD)
  print*, 'Non zero elements:', noneZero
  allocate(arrNZIndices(noneZero))
  arrNZIndices=arrayOfNonZeroElemIndices(typeSourceD,noneZero)
  print*, 'Non zero elem`s indices:'
  print '(20i5)', arrNZIndices

!$ Чтение динамических(!) матриц процессов (изменен s91mk/pritmd.f)
  initAbsorbMatrix   = 0.0d0
  initFissionMatrix  = 0.0d0
  leakageMatrix      = 0.0d0
  absMatrix  = 0.0d0
  fisMatrix  = 0.0d0
  leakMatrix = 0.0d0

  do step = 1, numMatr
      call rdind2(17,'LMAT',leakMatrix,2*numOfZones,0,step)
      leakageMatrix(:,numMatr-step+1) = leakMatrix(:)
      do zone = 1, numOfZones
          call rdind2(17,'FMAT',fisMatrix,2*numOfZones,zone,step)
          call rdind2(17,'AMAT',absMatrix,2*numOfZones,zone,step)
          initFissionMatrix(zone,:,numMatr-step+1) = fisMatrix(:)
          initAbsorbMatrix(zone,:,numMatr-step+1)  = absMatrix(:)
      enddo
  enddo

  deallocate(absMatrix, fisMatrix, leakMatrix)

  lnTime(1:numMatr-1) = 0.d0
  lnTime(numMatr) = maxTime

!$Формирование лог.интервалов для матриц процессов (2х-точность)
  lnTimeStep = dlog(maxTime/minTime)/(numMatr-2)
  lnTimeSum = 0.d0
  do step = numMatr-1, 2, -1
      lnTimeSum = lnTimeSum + lnTimeStep
      lnTime(step) = maxTime*dexp(-lnTimeSum)
  enddo
  print 4, lnTime(1:numMatr)
  print 7

  return

4 format(/' List of ln-times for fission matrices'/,(10es14.6e2))
6 format(10es14.6e2)
7 format(' ')

end



subroutine timeFinder

  use STK

  cnt = 1
  prprint = 0

  allocate(cntStep(numMatr), left(numMatr), center(numMatr),          &
           right(0:numMatr), linToLn(numMatr))

  cntStep = 0
  linToLn = 0
  left    = 0.0d0
  center  = 0.0d0
  right   = 0.0d0

  if (prprint == 1) print 6

  do step = 1, steps

    linTime = linTime + dt

    timeInIf     = 0.0d0

    do
      if (lnTime(cnt) < linTime .and. linTime >= lnTime(cnt+1)) then
        if (cnt > numMatr-1) exit
        left(cnt)   = 1.0d0
        right(cnt)  = 0.0d0
        center(cnt) = 0.0d0
        linToLn(cnt) = step
        if (prprint == 1)                                             &
          print 3,'++++',step,cnt,lnTime(cnt),linTime,lnTime(cnt+1),  &
                  left(cnt),center(cnt),right(cnt),cntStep(cnt)
        cnt = cnt + 1
      else
          exit
      endif
    enddo

    if (lnTime(cnt) < linTime .and. linTime < lnTime(cnt+1)) then
      left(cnt) = (linTime-lnTime(cnt))/(lnTime(cnt+1)-lnTime(cnt))
      if ((lnTime(cnt+1)-linTime) < dt) then
        center(cnt) = 0.0d0
        right(cnt)  = 1.0d0 - left(cnt)
        linToLn(cnt) = step
        if (prprint == 1)                                             &
          print 3,'1111',step,cnt,lnTime(cnt),linTime,lnTime(cnt+1),  &
                  left(cnt),center(cnt),right(cnt),cntStep(cnt)
      else
        do while (timeInIf < lnTime(cnt+1))
          timeInIf = linTime + cntStep(cnt)*dt
          cntStep(cnt) = cntStep(cnt) + 1
          if (timeInIf >= maxTime) exit
        enddo
        cntStep(cnt) = cntStep(cnt) - 2
        timeInIf  = timeInIf - dt
        right(cnt)  = (lnTime(cnt+1)-timeInIf)/                       &
                                        (lnTime(cnt+1)-lnTime(cnt))
        center(cnt) = dt/(lnTime(cnt+1)-lnTime(cnt))
        linToLn(cnt) = step
        if (prprint == 1)                                             &
          print 3,'    ',step,cnt,lnTime(cnt),linTime,lnTime(cnt+1),  &
                  left(cnt),center(cnt),right(cnt),cntStep(cnt)
      endif
      cnt = cnt + 1
    endif

  enddo

  if ((center(cnt-1) - right(cnt-1)) < 1.d-08) then
    cntStep(cnt-1) = cntStep(cnt-1) + 1
    if (prprint == 1) print*, 'last', cntStep(cnt-1), cnt-1
  endif

3 format(a4,2(i9,2x),2(es14.6e2,'  ::'),es14.6e2,'  |',3es14.6e2,i9)
4 format(10es14.6e2)
6 format(9x,' lin step ',' matrix ',2x,'left bound',9x,'lin step',    &
         9x,'right bound',8x,'left',10x,'center',8x,'right')

  return

end



!$    Solver for init fission source
subroutine initSource

  use STK

  real*8  :: sumR, sumLC
  integer :: startIndex, endIndex, m, mt

  allocate(matrStep(0:steps), signStep(steps),                        &
           normSource(numOfZones,steps))

  matrStep = 0
  signStep = 0
  prprint  = 0

  print*, ' '

  do step = 1, numMatr-1

    startIndex = linToLn(step)
    endIndex   = linToLn(step) + cntStep(step)

    do i = startIndex, endIndex
      if (cntStep(step) > 0) then
        if (i == startIndex) then
          signStep(i) = 1
          matrStep(i) = step

          if (prprint == 1) print 3, i, left(step), 0.0d0,            &
                    right(step), linToLn(step), cntStep(step), step
        elseif (i == endIndex) then
          signStep(i) = 2
          matrStep(i) = step

          if (prprint == 1) print 3, i, left(step), 0.0d0,            &
                    right(step), linToLn(step), cntStep(step), step
        else
          signStep(i) = 2
          matrStep(i) = step

          if (prprint == 1) print 3, i, 0.0d0, center(step), 0.0d0,   &
                                 linToLn(step), cntStep(step), step
        endif
      else
        signStep(i) = 3
        matrStep(i) = step

        if (prprint == 1) print 3, i, left(step), center(step),       &
                    right(step), linToLn(step), cntStep(step), step
      endif
    enddo

  enddo


  print*, ' '

  normSource = 0.0d0

  do i = 1, steps

    mt = matrStep(i)
    do j = 1, numOfZones

      sumLC = skipZeroSum(initFissionMatrix(:,j,mt),typeSourceD)
      sumR  = skipZeroSum(initFissionMatrix(:,j,mt-1),typeSourceD)

      if (signStep(i) == 1) then
        normSource(j,i) = left(mt)*sumLC + right(mt-1)*sumR
      endif

      if (signStep(i) == 2) then
        normSource(j,i) = center(mt)*sumLC
      endif

      if (signStep(i) == 3) then
        if (matrStep(i) - matrStep(i-1) > 1) then
          do m = matrStep(i-1) + 1, matrStep(i)
            sumLC = skipZeroSum(initFissionMatrix(:,j,m),typeSourceD)
            sumR  = skipZeroSum(initFissionMatrix(:,j,m-1),typeSourceD)
            normSource(j,i) = normSource(j,i) + sumLC*                &
                             (left(m) + center(m)) + right(m-1)*sumR
          enddo
        else
          normSource(j,i) = normSource(j,i) + (left(mt) +             &
                                  center(mt))*sumLC + right(mt-1)*sumR
        endif

      endif

    enddo

    if (prprint == 1) print 4, signStep(i), mt, i, sum(normSource(:,i))

  enddo

1 format(10es14.6e2)
4 format(3i8,10es14.6e2)
2 format(20i5)
3 format(i6,3es14.6e2,3i6)

  return

end



subroutine findSolve

  use STK

  real*8  :: fisRate, fisRateR, absRate, absRateR, leakRate, leakRateR
  real*8  :: Kcr, timeCounter
  integer :: m, mt, ind, n, stepForSumEnd, stepForSumStart, kk, j
  integer*8 :: k
  logical :: op

  allocate(fissionRate(numOfZones,0:steps), absorbRate(steps),        &
           leakageRate(steps))

  timeCounter = 0.0d0
  fisRate     = 0.0d0
  fisRateR    = 0.0d0
  absRate     = 0.0d0
  absRateR    = 0.0d0
  leakRate    = 0.0d0
  leakRateR   = 0.0d0

  fissionRate = 0.0d0
  absorbRate  = 0.0d0
  leakageRate = 0.0d0

  print 2
  print 3

  inquire(unit=chnl_out, opened=op)

!$ Main cycle per neutron generation
  do n = 1, neutronGenerationCount

  do i = 1, steps

    stepForSumStart = 1
    stepForSumEnd   = i
    if (n > 1) then
      stepForSumStart = 1 + i
      stepForSumEnd   = steps + i
      fissionRate(:,i) = 0.0d0
    endif


    do k = stepForSumStart, stepForSumEnd

      if (n == 1) ind = i - k + 1
      if (n > 1)  ind = steps + (i - k + 1)

      if (k <= steps) then
        mt = matrStep(k - 1)
        kk = k
      elseif (k > steps) then
        mt = matrStep(k - 1 - steps)
        kk = k - steps
      endif

!$    Fission & absorbtion & leak
      do j = 1, numOfZones

        fisRate =                                                     &
           skipZeroSum(initFissionMatrix(:,j,mt),fissionRate(:,ind))
        absRate =                                                     &
           skipZeroSum(initAbsorbMatrix(:,j,mt),fissionRate(:,ind))
        leakRate =leakageMatrix(j,mt)*fissionRate(j,ind)

        if (right(mt-1) /= 0.0d0 .and. mt > 0) then
          fisRateR =                                                  &
            skipZeroSum(initFissionMatrix(:,j,mt-1),fissionRate(:,ind))
          absRateR =                                                  &
            skipZeroSum(initAbsorbMatrix(:,j,mt-1),fissionRate(:,ind))
          leakRateR=leakageMatrix(j,mt-1)*fissionRate(j,ind)
        endif

        if (signStep(kk) == 1) then
          fissionRate(j,i) = fissionRate(j,i) + left(mt)*fisRate +    &
                                                 right(mt-1)*fisRateR
          absorbRate(i)  = absorbRate(i) + left(mt)*absRate +         &
                                                 right(mt-1)*absRateR
          leakageRate(i) = leakageRate(i) + left(mt)*leakRate +       &
                                                right(mt-1)*leakRateR
        endif

        if (signStep(kk) == 2) then
          fissionRate(j,i) = fissionRate(j,i) + center(mt)*fisRate
          absorbRate(i)  = absorbRate(i) + center(mt)*absRate
          leakageRate(i) = leakageRate(i) + center(mt)*leakRate
        endif

        if (signStep(kk) == 3) then
          if (matrStep(kk-1) - matrStep(kk-2) > 1) then
            do m = matrStep(kk-2) + 1, matrStep(kk-1)

              fisRate = skipZeroSum(initFissionMatrix(:,j,m),         &
                                                 fissionRate(:,ind))
              absRate = skipZeroSum(initAbsorbMatrix(:,j,m),          &
                                                 fissionRate(:,ind))
              leakRate =leakageMatrix(j,m)*fissionRate(j,ind)

              if (right(m-1) /= 0.0d0) then
                fisRateR = skipZeroSum(initFissionMatrix(:,j,m-1),    &
                                                 fissionRate(:,ind))
                absRateR = skipZeroSum(initAbsorbMatrix(:,j,m-1),     &
                                                 fissionRate(:,ind))
                leakRateR=leakageMatrix(j,m-1)*fissionRate(j,ind)
              endif

              fissionRate(j,i) = fissionRate(j,i) + (left(m) +        &
                          center(m))*fisRate + right(m-1)*fisRateR
              absorbRate(i) = absorbRate(i) + (left(m) +              &
                          center(m))*absRate + right(m-1)*absRateR
              leakageRate(i) = leakageRate(i) + (left(m) +            &
                          center(m))*leakRate + right(m-1)*leakRateR

            enddo

          else
            fissionRate(j,i) = fissionRate(j,i) + (left(mt) +         &
                        center(mt))*fisRate + right(mt-1)*fisRateR
            absorbRate(i) = absorbRate(i) + (left(mt) +               &
                        center(mt))*absRate + right(mt-1)*absRateR
            leakageRate(i) = leakageRate(i) + (left(mt) +             &
                        center(mt))*leakRate + right(mt-1)*leakRateR
          endif
        endif

      enddo

    enddo

    if (n == 1) fissionRate(:,i) = fissionRate(:,i) + normSource(:,i)

    fRate = sum(fissionRate(:,i))
    aRate = absorbRate(i)
    lRate = leakageRate(i)

    Kcr = fRate/(aRate + lRate)
    timeCounter = timeCounter + dt

    if (mod(i,modStep) == 0) then
      print 1, timeCounter, fRate, aRate, lRate, Kcr
      if (op) write(chnl_out, 1) timeCounter, fRate, aRate, lRate, Kcr

      if (timeCounter >= stopTime - 1.0d-12) then
        print*, 'Normal exit: according "STOP"'
        return
      endif

    endif

  enddo
  print 3

  absorbRate  = 0.0d0
  leakageRate = 0.0d0

  enddo


  inquire(unit=chnl_out, opened=op)
  if (op) close(chnl_out)


1 format(es14.6e2,'  |',3(es14.6e2,2x),' | ',es14.6e2)
2 format(5x,'time, s',7x,'fission rate',3x,'absorbtion rate',4x,      &
         'leak rate',9x,'K-critical')
3 format('------------------------------------------------------------&
  ----------------------')

  return
end
