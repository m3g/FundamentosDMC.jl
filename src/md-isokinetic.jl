"""

```
md_isokinetic(x::Vector{T},opt) where T
```

Performs a MD simulation starting from `x`, with options given by `opt`, with
an isokinetic bath. 

"""
function md_isokinetic(x0::Vector{T},opt::Options{T}) where T

  println("""
  -----------------------------------------------------------------------------------------
  NVT - Isokinetic-bath simulation 
  -----------------------------------------------------------------------------------------
  """)

  # aliases to simplify code
  n = length(x0)
  @unpack dt, nsteps, kavg_target, ibath, iequil = opt
  u(x) = potential(x,opt)
  f!(x,f) = forces!(x,f,opt)

  # Copy the initial point, to preserve it
  x = copy(x0)

  # Obtain initial velocities
  v = velocities(n,opt)

  # Open trajectory file for writting
  trajectory_file = open(opt.trajectory_file,"w")

  # log matrix will contain potential, kinetic, total energies, and "temperature" 
  log = zeros(nsteps,4)
 
  println(" Potential energy at initial point: " , u(x))
  println(" Kinetic energy at initial point: ", kinetic(v))
  println(" Total initial energy = ", u(x) + kinetic(v))

  # Write coordinates to trajectory file, and update log vectors
  printxyz(0.,x,opt,trajectory_file)
  log[1,:] .= (u(x), kinetic(v), u(x) + kinetic(v), kinetic(v)/n) 
  
  # Initialize velocity vector and save first set of forces 
  f = zeros(T,n)
  f!(f,x)
  flast = copy(f)

  # Running simulation
  time = 0.
  for istep = 1:nsteps

    # Write coordinates to trajectory file, and update log vectors

    # Updating positions 
    @. x = x + v*dt + 0.5*f*dt^2

    # Updating the forces
    @. flast = f
    f!(f,x)

    # Updating velocities
    @. v = v + 0.5*(f + flast)*dt 

    # Update step and print data
    time += dt
    kstep = kinetic(v) 
    ustep = u(x)
    energy = kstep + ustep 
    kavg = kstep / n

    # Stop if simulation exploded
    if ustep > 1e10
      println(" Simulation exploded: Energy = ", energy)
      close(trajectory_file)
      return log[1:istep,:]
    end

    # Save point to log
    log[istep,:] .= (ustep,kstep,energy,kavg)
    if mod(istep,opt.iprint) == 0 
      @printf(" TIME = %12.3f U = %12.5e K = %12.5e TOT = %12.5e \n", 
        time, ustep, kstep, energy
      )
      printxyz(time,x,opt,trajectory_file)
      log[istep,:] .= (ustep,kstep,energy,kavg)
    end

    #
    # Isokinetic-bath: rescale velocities to obtain target average kinetic energy
    #
    if istep <= iequil && mod(istep,ibath) == 0
      @. v = v * sqrt(kavg_target/kavg)
    end

  end

  close(trajectory_file)
  return log
end

