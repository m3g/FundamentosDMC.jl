
# Velocity Distribution

The velocities of particles in a (three-dimensional) molecular system follow the Maxwell-Boltzmann distribution:

$$f(v) ~dv= \left(\frac{m}{2 \pi kT}\right)^{3/2}\, 4\pi v^2 e^{ -\frac{mv^2}{2kT}} ~ dv$$

Essentially, this equation says that the probability that a particle has a velocity in the interval $v\pm dv$ is proportional to $v^2 e^{-v^2/kT}$. This is due to two reasons: first, the probability that a particle has energy $E$ is given by the Boltzmann distribution, $e^{-E/kT}$. Since the kinetic energy is $mv^2/2$, the square of the velocity appears in the exponent. Second, the *degeneracy* increases with $v^2$ as velocity increases, because the volume of space containing vectors of length $v$ increases with $v^2$ in three dimensions. This generates the pre-exponential term proportional to $v^2$.

Let us see if this actually occurs in our simulations.

## 10.1. Initial velocity distribution

We can plot the Maxwell-Boltzmann distribution to see its shape as a function of temperature:

```julia-repl
julia> mb(v,kT) = (1/(2π*kT))^(3/2)*4π*v^2*exp(-v^2/(2*kT))

julia> v = 0:0.1:5 # velocity range

julia> y = mb.(v,0.6) # the . applies to all elements of x

julia> plot(v,y,label="Maxwell-Bolztmann",xlabel="v",ylabel="frequency")
```

An interesting aspect to observe in simulations is how the velocity distribution evolves and
converges. In our programs
we have three different ways to generate the *initial* velocity distribution: with norms uniformly distributed between 0 and 1, with a normal distribution of norms, or all zero. This can be set with the `initial_velocities` parameter of `Options`. For example, initial velocities with different distributions can be generated as follows:

```julia-repl
julia> sys = System(n=100_000,sides=[50,50,50])
System{Point2D}:
 Number of particles = 100000
 Box sides = [100.0, 100.0]

julia> v_normal = init_velocities(sys,Options(initial_velocities=:normal))
10000-element Vector{Point2D}:
 [0.04326001565583558, 0.4435409371990999]
 ⋮
 [-0.13167046004485733, -0.8914869714395789]

julia> v_flat = init_velocities(sys,Options(initial_velocities=:flat))
10000-element Vector{Point2D}:
 [-0.6634658506319203, 0.9721022278467606]
 ⋮
 [1.773444140352424, 0.27389840410745514]
```

To see the velocity distribution in each case, we will use the `velocity_distribution` function, which reads the velocities stored in a trajectory file (the reason for this will become clear shortly):

```julia-repl
julia> printxyz(v_normal,sys,"vnormal.xyz")

julia> printxyz(v_flat,sys,"vflat.xyz")

julia> d_normal = velocity_distribution(sys,"vnormal.xyz");
Number of frames read = 1
 Average velocity = [-3.4000000016826124e-8, 4.899999999760096e-8]
 Average velocity norm = 0.9699959035972906
 Average kinetic energy = 0.5998727466282452

julia> d_flat = velocity_distribution(sys,"vflat.xyz");
Number of frames read = 1
 Average velocity = [-4.9000000006849124e-8, -1.2999999991358812e-8]
 Average velocity norm = 0.9506616041200406
 Average kinetic energy = 0.5999232099558015
```

Note that we generated `100_000` points to have a well-sampled distribution.
Both distributions have the same average kinetic energy, of course, corresponding to the desired temperature. The distributions, on the other hand, are very different:

```julia-repl
julia> plot(
           [d_normal,d_flat],
           label=["Normal" "Flat"],
           xlabel="velocity norm",
           ylabel="frequency",
           linewidth=2
       )
```
which results in:
```@raw html
<center>
<img src=../figures/velocities.svg>
</center>
```

The normal distribution is the one that resembles the Maxwell-Boltzmann distribution and is the default in the simulations.

## 10.2. Velocities at equilibrium

We will compare the initial distributions with a distribution obtained from
a simulation. For this, we will repeat the Langevin simulation, but now
in 3 dimensions, starting from zero velocities. Note that the simulation also saves the velocity trajectory:
```julia-repl
julia> sys = System(n=100,sides=[50,50,50])
System{Point3D}:
 Number of particles = 100
 Box sides = [50.0, 50.0, 50.0]

julia> minimize!(sys)

julia> md_out = md_langevin(
           sys,
           Options(
               lambda=0.01,
               nsteps=100_000,
               velocity_file="vel.xyz",
               initial_velocities=:zero
           )
       )
```

We compute the velocity distribution of this trajectory with:
```julia-repl
julia> v_sim = velocity_distribution(sys,"vel.xyz");
Number of frames read = 10001
 Average velocity = [0.003014992500749613, -0.04277096828317029, -0.01629199255074569]
 Average velocity norm = 1.2234527952231138
 Average kinetic energy = 0.8819808992851449
```

Plot this velocity distribution:
```julia-repl
julia> plot(
           v_sim,
           label="Langevin",
           xlabel="velocity norm",
           ylabel="frequency",
           linewidth=2
       )
```

Visually compare the result with what is expected for the Maxwell-Boltzmann distribution.

## 10.3. Fitting the data

The `v_sim` curve can be fitted to see how well, or poorly, it corresponds to the Maxwell-Boltzmann distribution. The general form is

$$f(v)dv = A v^2 e^{-v^2/2kT} dv$$

This function has only two adjustable parameters: $A$ and $kT$ (or simply $T$, but here we keep the constant together with the temperature).

We can fit this curve to the data using fitting software, or for example the [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) package:

### Fitting the model

```julia-repl
julia> using LsqFit

julia> @. f(v,p) = p[1]*v^2*exp(-v^2/(2*p[2]))

julia> x = v_sim[1]; y = v_sim[2]; # to make code clearer

julia> p0 = rand(2) # initial point

julia> fit = curve_fit(f, x, y, p0)
```

And we see how well the data are fitted by the model by computing the sum of squared residuals:

```julia-repl
julia> sum(fit.resid.^2)
2.016208460480086e-6
```

And we can visualize the fitted model graphically with:
```julia-repl
julia> pars = coef(fit)

julia> yfit = [ f(xi,pars) for xi in x ];

julia> plot(
           [ (x,y), (x,yfit) ],
           label=[ "Simulation"  "Fit" ],
           xlabel="velocity norm",
           ylabel="frequency",
           linewidth=2
       )
```

```@raw html
<center>
<img src=../figures/velocities_fit.svg>
</center>
```

Note that the error is very low, indicating that the simulation velocities fit the Maxwell-Boltzmann distribution well.
Furthermore, the second fit parameter should be the $kT$ that was defined for the simulation:
```julia-repl
julia> coef(fit)[2]
 0.5852713287386224
```

## 10.4. Complete summary code

```julia
using FundamentosDMC, Plots, LsqFit

# system setup
sys = System(n=100,sides=[50,50,50])
minimize!(sys)

# simulation
md_out = md_langevin(
    sys,
    Options(
        lambda=0.01,
        nsteps=100_000,
        velocity_file="vel.xyz",
        initial_velocities=:zero
    )
)

# plot velocity distribution
v_sim = velocity_distribution(sys,"vel.xyz");
plot(
    v_sim,
    label="Langevin",
    xlabel="velocity norm",
    ylabel="frequency",
    linewidth=2
)

# Fit model
@. f(v,p) = p[1]*v^2*exp(-v^2/(2*p[2]))
x = v_sim[1]; y = v_sim[2]; # to make code clearer
p0 = rand(2) # initial point
fit = curve_fit(f, x, y, p0)

# compare fit with simulation
pars = coef(fit)
yfit = [ f(xi,pars) for xi in x ]
plot(
    [ (x,y), (x,yfit) ],
    label=[ "Simulation"  "Fit" ],
    xlabel="velocity norm",
    ylabel="frequency",
    linewidth=2
)
savefig("./velocities_fit.pdf")

# Check kT
coef(fit)[2]
```
