!
! histogram: From an ascii data files, builds an histogram
!
! Usage: histogram data.dat [options]
! Options: -c [column to be considered]
!          -min [minimum value to be considered] 
!          -max [maximum value to be considered] 
!          -step [step of summation]
!          -nsteps [number of steps]
! 
!          Normalization (default: -int0):
!          -int0 : normalize to unity the integral of the histogram. 
!          -int1 : normalize the maximum probability to unity
!          -int2 : normalize such that each point contains the actual probablity
!                  of a value in the interval [ x-step, x+step ]
!
! All parameters except the input file name (data.dat) are optional. Default values are:
!
! column: 1
! minimum value: The mininum value from data.
! maximum value: The maximum value from data.
! step: (max - min) / nsteps
! number of steps: 100
!
! L. Martinez, IFSC-USP, Sep 6, 2011.
!
! Version 15.167
!

program histogram

implicit none
real :: min, max
integer :: col, nsteps, narg, iargc, ioerr, i, ndata, iline, first, idata
real :: val, step, average, sd, integral, pmax
real, allocatable :: xhist(:), yhist(:)
character(len=200) :: record, filename
character(len=17) :: number
logical int0, int1, int2

! Default parameter values

col = 1
min = +1.d30
max = -1.d30
nsteps = 100
step = -1.d0
first = 1
int0 = .false.
int1 = .false.
int2 = .false.

! Reading command line arguments

narg = iargc()
if ( narg == 0 ) then
  call print_options
  stop
end if

call getarg(1,record)
filename = trim(record)

i = 1
do while( i < narg )
  i = i + 1
  call getarg(i,record)
  select case (record)
    case('-c') 
      i = i + 1
      call getarg(i,record)
      read(record,*) col
    case('-max')
      i = i + 1
      call getarg(i,record)
      read(record,*) max
    case('-min')
      i = i + 1
      call getarg(i,record)
      read(record,*) min
    case('-nsteps')
      i = i + 1
      call getarg(i,record)
      read(record,*) nsteps
    case('-step')
      i = i + 1
      call getarg(i,record)
      read(record,*) step
    case('-first')
      i = i + 1
      call getarg(i,record)
      read(record,*) first
    case('-int0')
      int0 = .true.
    case('-int1')
      int1 = .true.
    case('-int2')
      int2 = .true.
    case default
      write(*,*) ' ERROR: Some option set incorrectly. '
      call print_options
      stop
  end select
end do
if ( (.not.int0).and.(.not.int1).and.(.not.int2) ) int0 = .true.

! If user defined maximum and minimum values to be used were not set, compute them
! from data file

open(10,file=filename,action='read')
ndata = 0
iline = 0
idata = 0
do
  read(10,"( a200 )",iostat=ioerr) record
  if ( ioerr /= 0 ) exit
  iline = iline + 1

  record = trim(record)
  if( record(1:1) == "#" ) cycle

  read(record,*,iostat=ioerr) (val, i = 1, col - 1), val
  if ( ioerr /= 0 ) cycle

  write(number,"( f17.10 )") val
  read(number,*,iostat=ioerr) val
  if ( ioerr /= 0 ) then
    write(*,"('# Warning: Value in line ', i8, ' is strange: ', e17.10, ': ignoring.')") &
             iline, val 
    cycle
  end if

  idata = idata + 1
  if ( idata >= first ) then
    ndata = ndata + 1
    min = dmin1(val,min)
    max = dmax1(val,max)
  end if

end do

if ( ndata == 0 ) then
  write(*,*) ' ERROR: Could not read any data (ndata = 0) '
  stop
end if

if ( step < 0.d0 ) step = (max - min) / nsteps

! Writting header and parameters 

write(*,"( '# Data file name: ', a )" ) trim(filename)
write(*,"( '# Number of data points: ', i10 )") ndata
write(*,"( '# Column: ', i10 )" ) col
write(*,"( '# Minimum x to be considered: ', f17.10 )" ) min
write(*,"( '# Maximum x to be considered: ', f17.10 )" ) max
write(*,"( '# Number of steps: ', i10 )" ) nsteps
write(*,"( '# Length of step: ', f17.10)") step

! Allocating the histogram arrays

allocate( xhist(nsteps+1), yhist(nsteps+1) )

! Assigning the value of the x coordinate of each histogram step 

do i = 1, nsteps + 1
  xhist(i) = min + (i-1)*(max - min)/nsteps
end do

! Computing the histogram

rewind(10)
average = 0.d0
idata = 0
do
  read(10,"( a200 )",iostat=ioerr) record
  if ( ioerr /= 0 ) exit

  record = trim(record)
  if( record(1:1) == "#" ) cycle

  read(record,*,iostat=ioerr) (val, i = 1, col - 1), val
  if ( ioerr /= 0 ) cycle
  write(number,"( f17.10 )") val
  read(number,*,iostat=ioerr) val
  if ( ioerr /= 0 ) cycle

  idata = idata + 1
  if ( idata >= first ) then
    average = average + val
    do i = 1, nsteps + 1
      if ( val >= ( xhist(i) - step ) .and. &
           val <= ( xhist(i) + step ) ) yhist(i) = yhist(i) + 1.d0
    end do 
  end if

end do
average = average / float(ndata)

! Computing the standard deviation

rewind(10)
sd = 0.d0
idata = 0
do
  read(10,"( a200 )",iostat=ioerr) record
  if ( ioerr /= 0 ) exit

  record = trim(record)
  if( record(1:1) == "#" ) cycle

  read(record,*,iostat=ioerr) (val, i = 1, col - 1), val
  if ( ioerr /= 0 ) cycle
  write(number,"( f17.10 )") val
  read(number,*,iostat=ioerr) val
  if ( ioerr /= 0 ) cycle

  idata = idata + 1
  if ( idata > first ) then
    sd = sd + (abs(val - average))**2
  end if

end do
sd = sqrt( sd / float(ndata - 1) ) 
write(*,"( '# Average: ',f17.10,' Standard deviation: ',f17.10 )" ) average, sd

! Normalizing the histogram 

integral = 0.d0
pmax = 0.d0
do i = 1, nsteps + 1
  yhist(i) = yhist(i) / ndata
  if ( int0 ) integral = integral + yhist(i)*step
  if ( int1 ) pmax = dmax1(pmax,yhist(i))
end do
if ( int0 ) then
  do i = 1, nsteps + 1
    yhist(i) = yhist(i) / integral
  end do
end if
if ( int1 ) then
  do i = 1, nsteps + 1 
    yhist(i) = yhist(i) / pmax
  end do
end if

! Printing histogram

do i = 1, nsteps + 1
  write(*,"( f17.10, tr2, f17.10 )") xhist(i), yhist(i)
end do

end program histogram

!
! Subroutine that prints the usage options
!

subroutine print_options

write(*,*) ' '
write(*,*) ' Run with: '
write(*,*) ' '
write(*,*) ' histogram data.dat -c 2 '
write(*,*) ' '
write(*,*) ' where 2 is the column of data.dat to be considered. '
write(*,*) ' '
write(*,*) ' Additional options: ' 
write(*,*) '   -step [real] : Increment of histogram sections. '
write(*,*) '   -min [real] : Minimum value of the histogram scale. ' 
write(*,*) '   -max [real] : Maximum value of the histogram scale. '
write(*,*) '   -nsteps [integer]: Number of steps. '
write(*,*) '   -first [integer]: First data point to be considered. '
write(*,*) ' '
write(*,*) '   Normalization: (default: int0) '
write(*,*) '   -int0 : normalize to unity the integral of the histogram. '
write(*,*) '   -int1 : normalize the maximum probability to unity. '
write(*,*) '   -int2 : normalize such that each point contains the actual probability '
write(*,*) '           of a value in the interval [ x-step, x+step ] '
write(*,*) ' '

end subroutine print_options













