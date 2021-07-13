"""

```
md_berendsen(sys::System{T},opt::Options) where T
```

Performs a MD simulation with a Berendsen bath. 

"""
function md_berendsen(sys::System{T},opt::Options=Options()) where T

  println("""
  -----------------------------------------------------------------------------------------
  NVT - Berendsen-bath simulation 
  -----------------------------------------------------------------------------------------
  """)

  # aliases to simplify code
  @unpack n, x0 = sys
  @unpack dt, nsteps, kT, iequil, tau = opt
  u(x) = potential(x,sys,opt)
  f!(f,x) = forces!(f,x,sys,opt)

  # Copy the initial point, to preserve it
  x = copy(x0)

  # Obtain initial velocities
  v = init_velocities(sys,opt)

  # Open trajectory file for writting
  trajectory_file = open(opt.trajectory_file,"w")
  velocity_file = open(opt.velocity_file,"w")

  # log matrix will contain potential, kinetic, total energies, and "temperature" 
  out = zeros(nsteps,4)
 
  println(" Potential energy at initial point: " , u(x))
  println(" Kinetic energy at initial point: ", kinetic(v))
  println(" Total initial energy = ", u(x) + kinetic(v))

  # Write coordinates to trajectory file, and update log vectors
  printxyz(0.,x,v,sys,trajectory_file,velocity_file)
  out[1,:] .= (u(x), kinetic(v), u(x) + kinetic(v), kinetic(v)/n) 
  
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
    remove_drift!(v)

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
      return out[1:istep,:]
    end

    # Save point to log
    out[istep,:] .= (ustep,kstep,energy,kavg)
    if mod(istep,opt.iprint) == 0 
      @printf(" TIME = %12.3f U = %12.5e K = %12.5e TOT = %12.5e \n", 
        time, ustep, kstep, energy
      )
      out[istep,:] .= (ustep,kstep,energy,kavg)
    end
    if mod(istep,opt.iprintxyz) == 0
      printxyz(time,x,v,sys,trajectory_file,velocity_file)
    end

    #
    # Berendsen-bath
    #
    if istep <= iequil
      lambda = sqrt( 1 + (dt/tau)*((dim(T)*kT/2)/kavg-1) )
      @. v = v * lambda
    end

  end

  close(trajectory_file)
  close(velocity_file)
  println(" Wrote trajectory file: ", opt.trajectory_file)
  println(" Wrote velocities file: ", opt.velocity_file)
  return out 
end

