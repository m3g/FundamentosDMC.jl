
"""

```
mc(sys::System{T},opt=Options=Options()) where T
```

Peforms a Monte-Carlos simulation.

"""
function mc(sys::System{T},opt::Options=Options()) where T

  println("""
  -----------------------------------------------------------------------------------------
  Monte-Carlo simulation
  -----------------------------------------------------------------------------------------
  """)

  # aliases to simplify code
  @unpack x0 = sys
  @unpack nsteps, kT, alpha = opt
  u(x) = potential(x,sys,opt)

  # Copy the initial point, to preserve it
  x = copy(x0)
  xtrial = copy(x0)

  # Open trajectory file for writting
  trajectory_file = open(opt.trajectory_file,"w")

  ucurrent = u(x)
  println(" Potential energy at initial point: " , ucurrent)

  # log vector will contain potential energy
  out = zeros(1)

  # Write coordinates to trajectory file, and update log vectors
  printxyz(0,x,sys,trajectory_file)
  out[1] = ucurrent 

  # Running simulation

  last_print = 0
  last_printxyz = 0
  naccepted = 0
  for istep in 1:nsteps

    # Creating new positions
    @. xtrial = x + alpha*randn(T)

    # Computing the potential energy of this new configuration
    utrial = u(xtrial)

    # Check if this conformation is accepted
    if (utrial < ucurrent) || 
       (exp( -(utrial-ucurrent) / kT) > rand())

      naccepted += 1
      ucurrent = utrial
      @. x = xtrial
      push!(out,ucurrent)

      if istep - last_print >= opt.iprint
        @printf(" STEP= %10i U = %12.5f AR = %12.5f\n",
                istep, ucurrent, naccepted/istep)
        last_print = istep
      end
      if istep - last_printxyz >= opt.iprintxyz
        printxyz(naccepted,x,sys,trajectory_file)
        last_printxyz = istep
      end

    end

  end
  println(" Final acceptance rate = ", naccepted/nsteps)

  close(trajectory_file)
  println(" Wrote trajectory file: ", opt.trajectory_file)
  return out
end
