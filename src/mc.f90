
function mc(pars=Parameters())

  ! Simulation parameters

  nsteps = 500000
  iprint = 5

  ! k*T (equal to average kinetic energy for 2D system)

  kt = 0.6

  ! Variance of the perturbation of the coordinates

  perturbation = 0.1d0

  ! Get parameters from the user

  write(*,"( 75('-'),/,' Monte-Carlo Simulation ',/,75('-') )")

  write(*,"( a,$ )") ' Enter the variance of the perturbations: ' 
  read(*,*) perturbation

  write(*,"( a,$ )") ' Enter the number of steps: '
  read(*,*) nsteps
 
  ! Output options

  printxyz = .true.
  printxyz = .true.

  ! Start from random positions

  call initial(n,x)

  ! Trajectory file

  if ( printxyz ) open(10,file='trajectory.xyz')

  ! Energy data file 

  open(20,file='energies.dat')

  ustep = potential(n,x)
  write(*,*) 'Potential energy at initial point: ', ustep
  call printx(n,x,0.d0)

  ! Running simulation

  naccepted = 0
  do istep = 1, nsteps

    ! Creating new positions

    do i = 1, n

      xvar = normal(0.d0,perturbation)
      xtrial(i,1) = x(i,1) + xvar
      xtrial(i,1) = image(xtrial(i,1))

      xvar = normal(0.d0,perturbation)
      xtrial(i,2) = x(i,2) + xvar
      xtrial(i,2) = image(xtrial(i,2))
      
    end do

    ! Computing the potential energy of this new configuration

    utrial = potential(n,xtrial)

    ! Check if this conformation is accepted

    call random_number(random)

    if ( utrial < ustep .or. &
         exp( -1.d0*(utrial-ustep) / kt ) > random &
       ) then

      naccepted = naccepted + 1

      ustep = utrial
      do i = 1, n
        x(i,1) = xtrial(i,1)
        x(i,2) = xtrial(i,2)
      end do

      if ( mod(naccepted,iprint) == 0 ) then
        write(*,"( ' STEP= ', i10, ' U = ', f12.5, ' AR = ', f12.5 )")  &
              istep, ustep, dble(naccepted)/istep
        write(20,*) istep, ustep
        if ( printxyz ) call printx(n,x,dble(istep))
      end if

    end if

  end do

  close(20)
  if ( printxyz ) close(10)

end program mc
