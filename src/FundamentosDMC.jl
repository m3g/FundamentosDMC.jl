module FundamentosDMC

using Printf: @sprintf, @printf
using StaticArrays: FieldVector, @SMatrix
using FastPow: @fastpow
using Statistics: mean

export System, Options
export Point2D, Point3D
export minimize!
export init_velocities
export printxyz, readxyz
export md, md_isokinetic, md_berendsen, md_langevin, mc
export radial_distribution, velocity_distribution
export simulate_gui

include("./Points.jl")
include("./System.jl")
include("./Options.jl")

include("./image.jl")

include("./potential.jl")
include("./kinetic.jl")
include("./remove_drift.jl")
include("./forces.jl")

include("./minimize.jl")
include("./printxyz.jl")

include("./init_velocities.jl")
include("./md-simple.jl")
include("./md-isokinetic.jl")
include("./md-berendsen.jl")
include("./md-langevin.jl")
include("./mc.jl")

include("./radial_distribution.jl")
include("./velocity_distribution.jl")

"""
    simulate_gui(; n=100, sides=(100.0, 100.0), kind=:md)

Opens a graphical interface to interactively run and visualize, in 2D, the
simulations described in the tutorial: microcanonical MD (`:md`), MD with an
isokinetic bath (`:md_isokinetic`), a Berendsen bath (`:md_berendsen`), a
Langevin bath (`:md_langevin`), and Monte Carlo (`:mc`).

The interface shows the motion of the particles together with live plots of
the potential, kinetic and total energies, and of the temperature (average
kinetic energy per particle), and exposes all the relevant `System` and
`Options` parameters. Changing any parameter or the simulation type
restarts the simulation with the new settings.

This function is implemented as a package extension and requires `GLMakie`
to be loaded:

```julia-repl
julia> using FundamentosDMC, GLMakie

julia> simulate_gui()
```

"""
function simulate_gui(; kwargs...)
    error("""

    `simulate_gui` requires the GLMakie package to be loaded. Install and load it with:

        import Pkg; Pkg.add("GLMakie")
        using GLMakie

    and call `simulate_gui()` again.

    """)
end

end

