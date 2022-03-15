module STK

  real*8, allocatable  :: lnTime(:),                                  &
                          normSource(:,:),                            &
                          fissionRate(:,:),                           &
                          absorbRate(:),                              &
                          leakageRate(:),                             &
                          initAbsorbMatrix(:,:,:),  absMatrix(:),     &
                          initFissionMatrix(:,:,:), fisMatrix(:),     &
                          leakageMatrix(:,:),       leakMatrix(:),    &
                          left(:), center(:), right(:),               &
                          typeSourceD(:)
  real, allocatable    :: typeSource(:)
  integer, allocatable :: cntStep(:), linToLn(:), matrStep(:),        &
                          signStep(:), arrNZIndices(:)

  real*8               :: minTime, stopTime, maxTime, dt, linTime,    &
                          lnTimeStep, lnTimeSum, timeInIf,            &
                          fRate, aRate, lRate
  integer              :: steps, numMatr, step, cnt,                  &
                          prprint, numOfZones, zone, noneZero,        &
                          neutronGenerationCount, modStep

  real*8               :: time(4)
  integer              :: params(4)

  integer        :: chnl_inp = 30,                                    &
                    chnl_out = 31
  character(100) :: nameinp = 'kin',                                  &
                    nameout = 'kin.out'


  contains

  integer function countNonZeroVectorElem(vector)
    real*8, intent(in) :: vector(numOfZones)
    integer            :: nonZero

    nonZero = 0
    do i = 1, numOfZones
      if (vector(i) > 0.0d0) nonZero = nonZero + 1
    enddo
    countNonZeroVectorElem = nonZero

    return
  end


  function arrayOfNonZeroElemIndices(vector,nonZero) result(indices)
    real*8, intent(in)   :: vector(numOfZones)
    integer,intent(in)   :: nonZero
    integer, allocatable :: indices(:)
    integer              :: counter

    allocate(indices(nonZero))
    indices = 0
    counter = 0

    do i = 1, numOfZones
      if (vector(i) > 0.0d0) then
        counter = counter + 1
        indices(counter) = i
      endif
    enddo

    return
  end


  real*8 function skipZeroSum(rowMatrix,vector)
    real*8, intent(in) :: vector(numOfZones), rowMatrix(numOfZones)
    real*8             :: vec

    vec = 0.0d0
!$  !Sum only non zero vector`s indices
    do j = 1, noneZero
      vec = vec + rowMatrix(arrNZIndices(j))*vector(arrNZIndices(j))
    enddo
    skipZeroSum = vec
    return
  end

end
