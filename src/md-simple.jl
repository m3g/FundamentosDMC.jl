"""

```
md(x::Vector{T},opt) where T
```

Performs a simple MD simulation starting from `x`, with options given by `opt` 

"""
function md(x::Vector{T},opt::Options{T}) where T

  println("""
  -----------------------------------------------------------------------------------------
  Simple MD
  -----------------------------------------------------------------------------------------
  """)

  # aliases to simplify code
  n = length(x)
  @unpack dt, nsteps, kavg_target = opt
  u(x) = potential(x,opt)
  f!(x,f) = forces!(x,f,opt)

  # Initialize velocities with proper average kinetic energy
  v = rand(Point,n)
  kavg = kinetic(v) / n
  @. v = sqrt(kavg_target/kavg)*v

  # Open trajectory files for writting
  trajectory_file = open(opt.trajectory_file,"w")
  energies_file = open(opt.energies_file,"w")
  velocities_file = open(opt.velocities_file,"w")
 
  println(" Potential energy at initial point: " , u(x))
  println(" Kinetic energy at initial point: ", kinetic(v))
  println(" Total initial energy = ", u(x) + kinetic(v))
  
  # Initialize velocity vector and save first set of forces 
  f = zeros(Point,n)
  f!(f,x)
  flast = copy(f)

  # Running simulation
  time = 0.
  printxyz(time,x,opt,trajectory_file)
  for istep = 1:nsteps

    # Updating positions 
    @. x = x + v*dt + 0.5*f*dt^2

    time += dt
    kstep = kinetic(v) 
    ustep = u(x)
    energy = kstep + ustep 
    kavg = kstep / n
    if mod(istep,opt.iprint) == 0 
      @printf(" TIME = %12.5f U = %12.5f K = %12.5f TOT = %12.5f \n", 
        time, ustep, kstep, energy
      )
      if opt.printvel
        for i = 1:n
          println(velocities_file,norm2(v[i]))
        end
      end
    end
    if mod(istep,opt.iprintxyz) == 0 
      printxyz(time,x,opt,trajectory_file)
    end

    # Stop if simulation exploded
    if ustep > 1e10
      println(" Simulation exploded: Energy = ", energy)
      close(trajectory_file)
      close(velocities_file)
      close(energies_file)
      return nothing
    end

    # Updating the forces
    @. flast = f
    f!(f,x)

    # Updating velocities
    @. v = v + 0.5*(f + flast)*dt 

  end
  close(trajectory_file)
  close(velocities_file)
  close(energies_file)

  return nothing
end

