"""

```
md(ts,tprod;pars)
```

Performs a simple MD simulation with time-step `ts` and parameters given by `pars`. 

"""
function md(pars=Parameters())
  println("
           -----------------------------------------------------------------------------------------
                                                    Simple MD
           -----------------------------------------------------------------------------------------
  ")

  # unpack some Parameters to simplify code
  @unpack n, ts, nsteps, kavg_target = pars 

  # closure some parameters to simplify code
  potential(x) = potential(x,pars)
  kinetic(x) = kinetic(x,pars)
  forces!(x,f) = forces!(x,f,pars)

  # Initial positions
  x = initial(n)

  # Initialize velocities with proper average kinetic energy
  v = rand(Point,n)
  kavg = kinetic(v) / n
  @. v = sqrt(kavg_target/kavg)*v

  # Open trajectory files for writting
  trajectory_file = open(pars.trajectory_file,"w")
  energies_file = open(pars.energies_file,"w")
  velocities_file = open(pars.velocities_file,"w")
 
  println(" Potential energy at initial point: " , potential(x))
  println(" Kinetic energy at initial point: ", kinetic(v))
  eini = potential(x) + kinetic(v)
  kini = kinetic(v)
  println(" Total initial energy = ', eini")
  printxyz(time,x,pars,trajectory_file)
  
  # Save first set of forces 
  forces!(x,f)
  flast = copy(f)

  # Running simulation

  time = 0.
  for istep = 1:nsteps

    # Updating positions 
    for i in 1:n
      x[i] = x[i] + v[i]*dt + 0.5*f[i]*dt^2
      x[i] = image(x[i]) # keep minimum image for vizualization
    end

    time += dt
    kstep = kinetic(v) 
    ustep = potential(x)
    energy = kstep + ustep 
    kavg = kstep / n

    if mod(istep,iprint) == 0 
      @printf(" TIME = %12.5f U = %12.5f K = %12.5f TOT = %12.5f ", time, ustep, kstep, energy)
      if pars.printvel
        for i = 1:n
          println(velocities_file,v[i][1]^2 + v[i][2]^2)
        end
      end
    end
    if mod(istep,iprintxyz) == 0 
      printxyz(time,x,pars,trajectory_file)
    end

    # Stop if simulation exploded
    if ustep > 1.e10
      println(" Simulation exploded: Energy = ", energy)
      close(trajectory_file)
      close(velocities_file)
      close(energies_file)
      return nothing
    end

    # Updating the forces
    for i in 1:n
      flast[i] = f[i]
    end
    forces!(x,f)

    # Updating velocities
    for i in 1:n
      v[i] = v[i] + 0.5*(f[i] + flast[i])*dt
    end

  end
  close(trajectory_file)
  close(velocities_file)
  close(energies_file)

  return nothing
end

