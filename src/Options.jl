"""

Structure containing the options and parameters of the simulation.

### Examples
```julia-repl
julia> opt = Options()

julia> opt = Options(nsteps=2000)

```

"""
@with_kw struct Options
    # Time
    dt::Float64 = 0.05
    nsteps::Int = 2000

    # System parameters
    eps::Float64 = 1.0
    sig::Float64 = 2.0

    # Temperature bath parameters
    initial_velocities::Symbol = :normal
    kT::Float64 = 0.6
    ibath::Int = 10
    iequil::Int = 10
    tau::Int = 10
    lambda::Float64 = 0.1
    alpha::Float64 = 0.1

    # Output and printing options
    printxyz::Bool = true
    iprintxyz::Int = 10
    iprint::Int = div(nsteps, 100)
    minprint::Bool = false
    trajectory_file::String = "traj.xyz"
    velocity_file::String = "vel.xyz"
end
