"""

Structure containing the options and parameters of the simulation: box `sides`, and `eps` and `sig` parameters of the Lennard-Jones interactions,
and some additional printing options. 

### Example
```julia-repl
julia> opt = Options()

julia> opt = Options(sides=200.,nsteps=2000)

```

"""
@with_kw_noshow struct Options{T}
  # Time
  dt::Float64 = 0.1
  nsteps::Int = 2000

  # System parameters
  sides::T = Point2D(100,100)
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

# Constructor that allows setting dimension of sides dynamically
function Options(;kargs...)
  isides = findfirst(isequal(:sides),keys(kargs))
  if isides === nothing
    sides = Point2D(100,100) 
    return Options{Point2D}(;sides=sides,kargs...)
  else
    sides = values(kargs)[isides]
    kargs = filter(x->x[1] != :sides, kargs)
    @show kargs
    if length(sides) == 2
      sides = Point2D(sides...)
      return Options{Point2D}(;sides=sides,kargs...)
    elseif length(sides) == 3
      sides = Point3D(sides...)
      return Options{Point3D}(;sides=sides,kargs...)
    else
      error("Dimension of sides must be 2 or 3.")
    end
  end
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


