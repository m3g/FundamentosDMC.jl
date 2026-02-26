
# Langevin Thermostat

The [md-langevin.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-langevin.jl) function implements Langevin temperature control.

The Langevin thermostat assumes that each real particle of the system is immersed in a fluid of much smaller particles, which slow it down through friction and, at the same time,
exert random forces on the real particles. The
frictional deceleration is proportional to the velocity, and the
random collisions cause instantaneous changes in the velocities of
the real particles. The trajectory of a real particle is thus
propagated by modifying the forces and velocities. The forces are
modified by introducing friction,

$$\vec{F}(t) = -\nabla V[\vec{x}(t)] - \lambda \vec{v}(t)$$

And so that the random collisions from the bath particles are
instantaneous (causing instantaneous changes in velocities), one
option is to introduce them as velocity changes:

$$v(t+\Delta t) = v(t) + a(t)\Delta t + \sqrt{2\lambda kT \Delta t}\delta(t)$$

Where $\delta(t)$ is a Gaussian random variable with zero mean and
standard deviation 1. The relationship between the friction coefficient, $\lambda$,
which decelerates the particles, and the intensity of the stochastic collisions,
$\sqrt{2\lambda kT}$, is a result of the fluctuation-dissipation theorem, which describes Brownian motion.

# 6.1. Parameter control and thermostatization

The friction coefficient, $\lambda$, controls the behavior of a Langevin dynamics. In this case, the program initializes velocities to zero, so that the effect of the thermostat is highlighted. We will print the trajectory more frequently so that the frictional effect at the start is noticeable:

For example,
```julia-repl
julia> out = md_langevin(
           sys,
           Options(lambda=0.1,dt=0.05,nsteps=20_000,iprintxyz=5,initial_velocities=:zero)
       );

```
Run the program with various parameters, in particular these:

| $\lambda$ | $\Delta t$ |
|:---------:|:----------:|
| 0.001 |  0.05 |
| 0.01  |  0.05 |
| 10.   |  0.05 |
| 10.   |  0.001|
|       |       |

After each run, observe the energy plots and the
trajectories. Discuss whether the temperature reached the desired value (kinetic
energy equal to 60 units), and whether the total energy is on average
constant. Observe the particle motion in each trajectory.
The consistency of the thermostat depends on the integration time step, in
particular for large couplings.

## 6.2. Complete summary code

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_langevin(
    sys,
    Options(lambda=0.1,dt=0.05,nsteps=20_000,iprintxyz=5,initial_velocities=:zero)
)
plot(
    out,ylim=[-100,100],
    label=["Potential" "Kinetic" "Total" "Temperature"],
    xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```
