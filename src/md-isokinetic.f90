program md_isokinetic

  use ff
  implicit none
  integer, parameter :: n = 100
  integer :: i, nsteps, iprint, istep, ibath, nequil, iprintxyz, seed
  double precision :: x(n,2), v(n,2), f(n,2), time, dt, &
                      potential, kinetic, eini, energy, kini, kstep, &
                      ustep, tprint, tprod, tbath, tequil, flast(n,2), totaltime, &
                      kavg, kavg_target, tprintxyz, image, vx, vy, random
  logical :: restart, printvel
 
  ! Initialize random number generator

  call seed_from_time(seed)
  call init_random_number(seed)

  ! Simulation parameters

  time = 0.d0
  tequil = 500.d0
  tprod= 2000.d0
  tprint = 1.d0
  tprintxyz = 2.d0
  tbath = 50.d0
  printvel = .false.
  restart = .false.

  ! Get tbath from the user

  write(*,"( 75('-'),/,' MD with Isokinetic bath ',/,75('-') )")

  write(*,"( a,$ )") ' Enter the time between each velocity rescaling: ' 
  read(*,*) tbath

  write(*,"( a,$ )") ' Enter the equilibration time: '
  read(*,*) tequil
 
  ! Target average kinetic energy (equals to kT for 2D)

  kavg_target = 0.6

  ! Time-step

  dt = 0.05d0

  ! Initialize indices and coordinates

  totaltime = tequil + tprod
  nsteps = int( totaltime / dt )
  iprint = int( tprint / dt )
  iprintxyz = int( tprintxyz / dt ) 
  ibath = int( tbath / dt )
  nequil = int( tequil / dt )

  ! Initial point

  call initial(n,x)

  ! Initialize velocities are very small, to see thermoalization

  do i = 1, n
    call random_number(random)
    v(i,1) = -1.d0 + 2.d0*random
    call random_number(random)
    v(i,2) = -1.d0 + 2.d0*random
  end do
  kstep = kinetic(n,v)
  kavg = kstep / dble(n)
  do i = 1, n
    v(i,1) = v(i,1)*1.d-5*dsqrt(kavg_target/kavg)
    v(i,2) = v(i,2)*1.d-5*dsqrt(kavg_target/kavg)
  end do

  ! Trajectory file

  open(10,file='trajectory.xyz')

  ! Energy data file 

  open(20,file='energies.dat')

  ! Velocity distribution file

  open(30,file='velocities.dat')
 
  write(*,*) ' Potential energy at initial point: ', potential(n,x)
  write(*,*) ' Kinetic energy at initial point: ', kinetic(n,v)
  eini = potential(n,x) + kinetic(n,v)
  kini = kinetic(n,v)
  write(*,*) ' Total initial energy = ', eini
  call printx(n,x,time)
  call force(n,x,f)
  do i = 1, n
    flast(i,1) = f(i,1)
    flast(i,2) = f(i,2)
  end do

  ! Running simulation

  time = 0.d0
  do istep = 1, nsteps

    ! Updating positions 

    do i = 1, n
      x(i,1) = x(i,1) + v(i,1)*dt + 0.5d0*f(i,1)*(dt**2)
      x(i,2) = x(i,2) + v(i,2)*dt + 0.5d0*f(i,2)*(dt**2)
      x(i,1) = image(x(i,1))
      x(i,2) = image(x(i,2))
    end do

    time = time + dt
    kstep = kinetic(n,v) 
    ustep = potential(n,x)
    energy = kstep + ustep 
    kavg = kstep / dble(n)

    if ( mod(istep,iprint) == 0 ) then
      write(*,"( ' TIME= ', f12.5, ' U = ', f12.5, ' K = ', f12.5, ' TOT = ', f12.5 )") &
              time - tequil, ustep, kstep, energy
      write(20,*) time - tequil, ustep, kstep, energy
      if ( printvel ) then
        do i = 1, n
          write(30,*) v(i,1)**2 + v(i,2)**2
        end do
      end if
    end if
    if ( mod(istep,iprintxyz) == 0 ) call printx(n,x,time)

    ! Stop if simulation exploded

    if ( ustep > 1.d10 ) then
      write(*,*) ' Simulation exploded: Energy = ', energy
      stop
    end if

    ! Isokinetic bath

    if ( time < tequil .and. ( mod(istep,ibath) == 0 ) ) then

      ! Rescale velocities to target average kinetic energy

      vx = 0.d0
      vy = 0.d0
      do i = 1, n
        v(i,1) = v(i,1)*dsqrt(kavg_target/kavg)
        v(i,2) = v(i,2)*dsqrt(kavg_target/kavg)
        vx = vx + v(i,1)
        vy = vy + v(i,2)
      end do

      ! Remove drift
       
      vx = vx / n
      vy = vy / n
      do i = 1, n
        v(i,1) = v(i,1) - vx
        v(i,2) = v(i,2) - vy
      end do

    end if

    ! Updating the force

    do i = 1, n
      flast(i,1) = f(i,1)
      flast(i,2) = f(i,2)
    end do
    call force(n,x,f)

    ! Updating velocities
      
    do i = 1, n
      v(i,1) = v(i,1) + 0.5d0*(f(i,1)+flast(i,1))*dt
      v(i,2) = v(i,2) + 0.5d0*(f(i,2)+flast(i,2))*dt
    end do

  end do

  close(30)
  close(20)

end program md_isokinetic
