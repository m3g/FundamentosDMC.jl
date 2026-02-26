
# The Simulated System

The simulation consists of a fluid of 100 (monatomic) particles that
interact via a
Lennard-Jones potential, in a two-dimensional, periodic system.

$V = 4\varepsilon \left( \frac{\sigma^{12}}{r^{12}} - \frac{\sigma^6}{r^6} \right)$

Open the file [potential.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/potential.jl) and understand the implementation of the
potential energy calculation. Note that the calculation depends on 3
parameters: $\varepsilon$, $\sigma$, and the size of the periodic system. The
parameters are defined in the `opt` data structure, which is an input argument (we will see later how to use it).

The file [forces.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/forces.jl) contains the force calculation (the gradient
of the potential), and the file [kinetic.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/kinetic.jl) contains the
kinetic energy calculation. Since the system uses periodic boundary
conditions, coordinates must always be computed relative
to the minimum image. The minimum image calculation is implemented in
the file [image.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/image.jl). It is worthwhile to understand the
implementation of each of these functions, which are common to all
the methods we will describe.

## 2.1. Initial coordinates

To begin, open a `Julia` session and run:
```julia-repl
julia> using FundamentosDMC
```
The initial system can be created randomly using:
```julia-repl
julia> sys = System(n=100,sides=[100,100])
System{Point2D}:
 Number of particles = 100
 Box sides = [100.0, 100.0]

```
which generates `100` random points in 2 dimensions with coordinates between `[-50,-50]` and `[50,50]` in this case. `Point2D` is a variable type that represents a point in two dimensions. Later we will see that all the code is generic, and we can perform simulations in 3 dimensions simply by changing the type of variable used.

## 2.1. Simulation parameters and options

Simulation parameters are controlled via the `Options` structure. For example, to adjust the time step, we pass the `dt` parameter to the structure. This can be done in the call to the simulation functions, as we will see.

```julia-repl
julia> Options(dt=0.1)
Options
  dt: Float64 0.1
  nsteps: Int64 2000
  eps: Float64 1.0
  sig: Float64 2.0
  initial_velocities: Symbol normal
  kT: Float64 0.6
  ibath: Int64 10
  iequil: Int64 10
  tau: Int64 10
  lambda: Float64 0.1
  alpha: Float64 0.1
  printxyz: Bool true
  iprintxyz: Int64 10
  iprint: Int64 20
  minprint: Bool false
  trajectory_file: String "traj.xyz"
  velocity_file: String "vel.xyz"
```

In this case, we adjusted the time step manually and kept all other options at their default values. Each of these parameters will be discussed in due course. Note that they define the size, force field ($\varepsilon$ and $\sigma$), $kT$ (temperature), and the names of the output files.

## 2.3. Energy minimization

Next, the initial configuration will be modified using
the Gradient Descent method to minimize the energy. The
method consists of moving the particles according to the first-order Taylor
approximation, in the direction of decreasing energy:

$\vec{x}_{i+1} = \vec{x}_i - \nabla V(\vec{x}_i) \Delta x$

If the energy at point $\vec{x}_{i+1}$ is less than the energy at
point $\vec{x}_i$, the point $\vec{x}_{i+1}$ is accepted and the process is
repeated. If not, $\Delta x$ is reduced ($\Delta x = \Delta x / 2$), and
a new point $\vec{x}_{i+1}$ is computed. Since the approximation should be
a good approximation in the vicinity of the current point ($\vec{x}_i$), a
negative gradient guarantees that the function decreases for sufficiently
small $\Delta x$. The process is stopped when the norm of the
gradient is small, or when too many points have been tested. In mechanics, $-\nabla V = \vec{F}$, so the function that computes the gradient is the same one that computes the forces in the simulation. Open the file [minimize.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/minimize.jl) to
discuss how the initial configuration is created.

Before running the energy minimization, let us copy the initial configuration for comparison:
```julia-repl
julia> x0 = copy(sys.x0)
```

Then, we minimize the energy with the `minimize!` function:
```julia-repl
julia> minimize!(sys)
Energy before minimization: 38322.72337856496
Energy after minimization: -74.15646912098042
```

!!! note
    In `Julia`, it is a convention that functions which modify their arguments end with `!`. In this case, the function will modify the particle positions, `sys.x0`.

If you want to see the progression of energies, enable the `minprint` option with `minimize!(sys,Options(minprint=true))`.

We can quickly see what happened to the particles by plotting them. First, we generate a plot of the points before minimization:
```julia-repl
julia> using Plots

julia> scatter(Tuple.(x0))
```
The points should be randomly distributed, and in particular some points should be very close to others, generating very repulsive potentials.

Then, we plot the energy-minimized configuration:
```julia-repl
julia> scatter(Tuple.(sys.x0))
```
and you will notice that the points now have a new arrangement: there are points forming clusters, because the Lennard-Jones potential is attractive at long distances. But there are no longer very close points generating very large repulsions.

This energy-minimized initial configuration will be used in our simulations.

!!! note
    Plots can be saved as image files with `savefig("plot.pdf")`, for example.

## 2.4. Temperature

The system temperature is a parameter also defined internally
in the program (it can be modified at will, but we will not do so). The
temperature is defined from the average kinetic energy associated with each
degree of freedom of the system. In the case where all
motion can be written as translations, the definition is

$\frac{1}{2}kT = \left< \frac{1}{2} m v_x^2\right>$

where the average over $v_x$ here is equivalent to an average over
any other translational degree of freedom. In a three-dimensional system, therefore,

$\left<\frac{1}{2}m |\vec{v}|^2 \right> =
\left<\frac{1}{2}m \left(v_x^2 + v_y^2 + v_z^2\right) \right> =
3\left< \frac{1}{2} m v_x^2 \right> = \frac{3}{2}kT$

which is the usual result.

Our simulations are of a two-dimensional system. In this case,

$\left< \frac{1}{2}m |\vec{v}|^2 \right> =
\left< \frac{1}{2}m \left(v_x^2 + v_y^2\right)\right> =
2\left< \frac{1}{2}m v_x^2 \right> = kT$

In all the code it has been chosen to target simulating the system at
the temperature corresponding to $kT = 0.6$ units. The simulated systems have 100 particles and are two-dimensional, so the average kinetic energy is $100kT=60$ units. The initial velocities will be generated randomly at the start of the simulation.

## 2.5. Complete summary code

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
x0 = copy(sys.x0)
minimize!(sys)
scatter(Tuple.(x0))
scatter(Tuple.(sys.x0))
```

Note: the `Tuple.(x)` command converts the vector of vectors into a vector of pairs (tuples), which is correctly interpreted by `Plots` as a single data series.



