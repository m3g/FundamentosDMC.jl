"""

Structure containing the options and parameters of the simulation: box `sides`, and `eps` and `sig` parameters of the Lennard-Jones interactions,
and some additional printing options. 

### Example
```julia-repl
julia> opt = Options()

julia> opt = Options(side=200.,nsteps=2000)

```

"""
@with_kw struct Options
  # Time
  dt::Float64 = 0.1
  nsteps::Int = 2000

  # System parameters
  sides::Point = Point(100,100)
  eps::Float64 = 1.
  sig::Float64 = 2.

  # Temperature bath parameters
  kavg_target = 0.6
  ibath::Int = 1

  # Output and printing options
  printxyz::Bool = true
  printvel::Bool = false
  iprint::Int = 1
  iprintxyz::Int = 2
  trajectory_file::String = "traj.xyz"
  energies_file::String = "energies.dat"
  velocities_file::String = "velocities.dat"
end


