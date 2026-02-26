
# Berendsen Temperature Control

The [md-berendsen.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-berendsen.jl) function implements Berendsen temperature
control. This method is also based on velocity
scaling, but is smoother. Velocities are scaled by

$$\lambda = \left[
1 + \frac{\Delta t}{\tau} \left(
\frac{T_0}{T(t)} -1
\right)
\right]^{1/2}$$

where $\Delta t$ is the integration time step and $\tau$ is a parameter that
defines the rate at which scaling is applied. The
scaling is smoother and slower.

# 5.1. Parameter control and thermostatization

The parameter $\tau$ is adjusted with the `tau` option of `Options`. For example:

```julia-repl
julia> sys = System(n=100,sides=[100,100])

julia> minimize!(sys)

julia> out = md_berendsen(sys,Options(tau=50,iequil=500,nsteps=20_000));
```

Try different parameters, with `20_000` simulation steps. Among others, try these:

| $\tau$ |  $i_{\mathrm{equil}}$ |
|:------:|:---------------------:|
|   50   |  500 |
|   50   |  1500|
|  300   |  1500|
|  300   |  3000|
|        |      |

Observe the resulting energy plots, using the same commands as before:
```julia-repl
julia> using Plots

julia> plot(
           out,ylim=[-100,100],
           label=["Potential" "Kinetic" "Total" "Temperature"],
           xlabel="step"
       )
```

Observe the smoothness, or lack thereof, of the total energy curve. Check whether
the kinetic energy approached the desired average energy ($kT=60$).

## 5.2. Complete summary code

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_berendsen(sys,Options(tau=50,iequil=500,nsteps=20_000))
plot(
    out,ylim=[-100,100],
    label=["Potential" "Kinetic" "Total" "Temperature"],
    xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```

