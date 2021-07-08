"""

Structure containing the options and parameters of the simulation: box `sides`, and `eps` and `sig` parameters of the Lennard-Jones interactions,
and some additional printing options. 

### Example
```julia-repl
julia> opt = Options()

julia> opt = Options(sides=200.,nsteps=2000)

```

"""
struct Options{T}
  # Time
  dt::Float64
  nsteps::Int 

  # System parameters
  sides::T
  eps::Float64
  sig::Float64

  # Temperature bath parameters
  kavg_target
  ibath::Int
  iequil::Int

  # Output and printing options
  printxyz::Bool
  printvel::Bool
  iprint::Int
  iprintxyz::Int
  trajectory_file::String
end

function Options(;
  dt::Float64 = 0.05,
  nsteps::Int = 2000,
  sides = nothing,
  eps::Float64 = 1.,
  sig::Float64 = 2.,
  kavg_target = 0.6,
  ibath::Int = 1,
  iequil::Int = div(nsteps,10),
  printxyz::Bool = true,
  printvel::Bool = false,
  iprint::Int = 1,
  iprintxyz::Int = 2,
  trajectory_file::String = "traj.xyz",
)

  if sides === nothing
    sides = Point2D(100,100)
  else
    if length(sides) == 2
      sides = Point2D(sides...)
    elseif length(sides) == 3
      sides = Point3D(sides...)
    else
      error("Dimension of sides must be 2 or 3.")
    end
  end
  Options{typeof(sides)}(dt,nsteps,sides,eps,sig,kavg_target,
          ibath,iequil,printxyz,printvel,iprint,iprintxyz,
          trajectory_file)

end

function Base.show(io::IO,opt::Options)
  println(typeof(opt))
  println("-------------------")
  println("Simulation options:")
  print("-------------------")
  for (i,field) in pairs(fieldnames(Options))
    print("\n",fieldname(Options,i)," = ",getfield(opt,field))
  end
end


