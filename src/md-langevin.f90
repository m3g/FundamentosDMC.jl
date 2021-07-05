program md_langevin

  use ff
  implicit none
  integer, parameter :: n = 100
  integer :: i, nsteps, iprint, istep, iprintxyz, seed
  double precision :: x(n,2), v(n,2), f(n,2), flast(n,2), time, dt, &
                      potential, kinetic, eini, energy, kini, kstep, &
                      ustep, tprint, tprod, totaltime, &
                      kavg, kavg_target, tprintxyz, normal, image, lambda
  logical :: printvel
  character(len=200) :: record
 
  ! Initialize random number generator

  call seed_from_time(seed)
  call init_random_number(seed)

  ! Simulation parameters

  tprod = 2000.d0
  tprint = 1.d0
  tprintxyz = 2.d0
  printvel = .false.

  ! Get total time from command line if given

  if ( iargc() > 0 ) then
    call getarg(1,record)
    read(record,*) tprod
  end if

  ! Time-step

  dt = 0.05d0

  ! Friction parameter for Langevin bath

  lambda = 0.5d0

  ! Get parameters of the bath from the user

  write(*,"( 75('-'),/,' MD with Langevin bath ',/,75('-') )")

  write(*,"( a,$ )") ' Enter bath coupling constant (lambda): ' 
  read(*,*) lambda

  write(*,"( a,$ )") ' Enter time-step: '
  read(*,*) dt

  ! Target average kinetic energy (equals to kT for 2D)

  kavg_target = 0.6

  ! Initialize indices and coordinates

  totaltime = tprod
  nsteps = int( totaltime / dt )
  iprint = int( tprint / dt )
  iprintxyz = int( tprintxyz / dt ) 

  ! Initial coordinates

  call initial(n,x)

  ! Trajectory file

  open(10,file='trajectory.xyz')
  call printx(n,x,time)

  ! Here we will initialize with null velocities, to see thermalization

  do i = 1, n
    v(i,1) = 0.d0
    v(i,2) = 0.d0
  end do

  ! Energy data file 

  open(20,file='energies.dat')

  write(*,*) ' Potential energy at initial point: ', potential(n,x)
  write(*,*) ' Kinetic energy at initial point: ', kinetic(n,v)
  eini = potential(n,x) + kinetic(n,v)
  kini = kinetic(n,v)
  write(*,*) ' Total initial energy = ', eini
  call force(n,x,f)

  ! Running simulation

  time = 0.d0
  do istep = 1, nsteps

    ! Updating positions

    do i = 1, n
      x(i,1) = x(i,1) + v(i,1)*dt + 0.5d0*f(i,1)*dt**2 
      x(i,2) = x(i,2) + v(i,2)*dt + 0.5d0*f(i,2)*dt**2
      x(i,1) = image(x(i,1))
      x(i,2) = image(x(i,2))
    end do

    ! Updating deterministic forces

    do i = 1, n
      flast(i,1) = f(i,1)
      flast(i,2) = f(i,2)
    end do
    call force(n,x,f)

    ! Add friction force

    do i = 1, n
      f(i,1) = f(i,1) - lambda*v(i,1)
      f(i,2) = f(i,2) - lambda*v(i,2)
    end do

    ! Updating velocities; v = v + 0.5*(f1/m+f2/m)dt - lambda v dt + sqrt(2 lambda kT dt) * normal

    do i = 1, n
      v(i,1) = v(i,1) + 0.5d0*(f(i,1)+flast(i,1))*dt &
               + dsqrt(2.d0*lambda*kavg_target*dt)*normal(0.d0,1.d0)
      v(i,2) = v(i,2) + 0.5d0*(f(i,2)+flast(i,2))*dt &
               + dsqrt(2.d0*lambda*kavg_target*dt)*normal(0.d0,1.d0)
    end do

    time = time + dt
    kstep = kinetic(n,v) 
    ustep = potential(n,x)
    energy = kstep + ustep 
    kavg = kstep / dble(n)

    if ( mod(istep,iprint) == 0 ) then
      write(*,"( ' TIME= ', f12.5, ' U = ', f12.5, ' K = ', f12.5, ' TOT = ', f12.5 )") &
              time, ustep, kstep, energy
      write(20,*) time, ustep, kstep, energy
    end if
    if ( mod(istep,iprintxyz) == 0 ) call printx(n,x,time)

    ! Stop if simulation exploded

    if ( ustep > 1.d10 ) then
      write(*,*) ' Simulation exploded: Energy = ', energy
      stop
    end if

  end do

  close(20)
  close(10)

end program md_langevin

