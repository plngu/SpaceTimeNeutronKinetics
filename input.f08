subroutine input

!$*******************************************************************
!$ Чтение входного файла кинетики                                   *
!$                                                                  *
!$ Значение ключей:                                                 *
!$ MATRIX   - кол-во матриц делений (процессов): `int`              *
!$ INTERVAL - временной интервал решения уравнений: `float` `float` *
!$            от мин. до макс. времени жизни нейтрона               *
!$ LIN_STEP - шаг интегрирования: `float`                           *
!$ STOP     - завершить расчет на времени t: `float` &              *
!$ ZONES    - кол-во зон в XT26: `int`                              *
!$ PRINT    - кол-во печатей, mod(step, PRINT) == 0: `int` &        *
!$ OUTFILE  - включить запись во входной файл: `string` &           *
!$ N_GNR    - кол-во поколений нейтронов: `int`                     *
!$ END      - конец чтения                                          *
!$                                                                  *
!$ & - пометка необязательного ключа                                *
!$*******************************************************************

  use STK

  integer        :: ind, er
  character(8)   :: key
  character(100) :: str
  logical        :: ex

  time(1:4)  = [0.0d0, 0.0d0, 0.0d0, 1.0d16]
  params(1:4) = [0, 0, 0, 0]

  inquire (file=trim(nameinp), exist=ex)
  if (ex.eqv..false.) stop '# Error: input file "kin" doesn`t exist'
  open(chnl_inp, file=trim(nameinp), action='read')
  print'(/,a)','# Start reading input file of kinetic module'

  do

    ind = ind + 1

    read(chnl_inp, fmt=4, end=2) key, str
    print 4, key, str

    if ((key(1:1).eq.'!').or.(key(1:1).eq.'*')) cycle

    select case(trim(key))
      case('', 'BEGIN')
        cycle
      case('MATRIX')
        read(str,*,iostat=er) params(1)
      case('INTERVAL')
        read(str,*,iostat=er) time(1:2)
      case('LIN_STEP')
        read(str,*,iostat=er) time(3)
      case('STOP')
        read(str,*,iostat=er) time(4)
      case('ZONES')
        read(str,*,iostat=er) params(2)
      case('PRINT')
        read(str,*,iostat=er) params(3)
      case('OUTFILE')
        read(str,*,iostat=er) nameout
        open(chnl_out, file=trim(nameout), action='write')
      case('N_GNR')
        read(str,*,iostat=er) params(4)
      case('END')
        exit
      case default
        print 1, trim(key), ind
        stop
    end select

  enddo

!$ Check
  if (any(params <= 0)) stop 'Wrong keys values'
  do i = 1, size(time)
    if (time(i) <= 0.0d0) stop 'Wrong times values'
  enddo

1 format(' # Error: Unknown key ','"',a,'"',' in ',i5, ' string')
2 close(chnl_inp)
4 format(a8,a100)

  return
end
