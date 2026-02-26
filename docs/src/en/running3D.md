# Running in 3D

All the code here can generate three-dimensional simulations. In fact, the only code that has any
defined dimensionality is that which writes the output file (`printxyz`) and that which computes the radial distribution
function (`radial_distribution`). All other code is generic, and what it does depends only on the type of the input variable, i.e., on the dimension of the box and the points.

## 9.1. Initializing the system in 3 dimensions:

Simply create the system in a three-dimensional box:

```julia-repl
julia> sys = System(n=100,sides=[50,50,50])
System{Point3D}:
 Number of particles = 100
 Box sides = [50.0, 50.0, 50.0]
```
Note that the points are now of type `Point3D`. Indeed, the coordinates, `sys.x0`, are three-dimensional:

```julia-repl
julia> sys.x0
100-element Vector{Point3D}:
 [1.7209592667160578, 2.6345231665739135, 1.9123751433362468]
 [0.8615941909462954, -5.3152138387533565, -1.6775684653805492]
 [5.492878544337906, -7.041239125665051, -4.395671467880198]
⋮
 [-2.229487380572972, -2.4612177799986634, 8.742176399685668]
 [3.9099597073208123, 7.797739759305074, 4.3056981321251655]
```

The energy can then be minimized as before:
```julia-repl
julia> minimize!(sys)
Energy before minimization: 1.8502746330219281e6
Energy after minimization: -451.7109718153855
```

## 9.2. Running and visualizing the simulation

The simulation, whatever type it is, can be run normally. The exact same code is used. For example, a Langevin dynamics simulation can be run with:

```julia-repl
julia> md_out = md_langevin(
          sys,
          Options(lambda=0.01,nsteps=20_000,trajectory_file="md.xyz")
       )
```

The simulation can be visualized in `VMD` with:
```
shell> vmd md.xyz
```

And we can visualize the radial distribution function with:

```julia-repl
julia> rmd, gmd = radial_distribution(sys,"md.xyz")

julia> plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")

```

## 9.3. Complete summary code

```julia
using FundamentosDMC, Plots

sys = System(n=100,sides=[50,50,50])
minimize!(sys)

md_out = md_langevin(
    sys,
    Options(lambda=0.01,nsteps=20_000,trajectory_file="md.xyz")
)
rmd, gmd = radial_distribution(sys,"md.xyz")
plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
```

