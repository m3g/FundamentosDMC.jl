
# Isokinetic Temperature Control

The function implemented in [md-isokinetic.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-isokinetic.jl) implements isokinetic temperature control. In this method, velocities are scaled by a
parameter $\lambda = \sqrt{T_0/T}$ at regular intervals, to
thermostatize the system at temperature $T_0$.
Its execution requires the definition of two additional
parameters: the time interval between two velocity scalings,
and the equilibration time. The equilibration time is
the time during which the scalings are performed. The
goal should be to obtain a stable simulation, with average kinetic energy
appropriate to the desired value (60 units), after equilibration.

# 4.1. Parameter control and thermostatization

The system is initialized in the same way as before:
```julia-repl
julia> using FundamentosDMC, Plots

julia> sys = System(n=100,sides=[100,100])

julia> minimize!(sys);
```

Then, we run the simulation, now with an isokinetic thermostat, for `2000` steps, of which `200` are equilibration steps. The thermostat is applied every `ibath` steps:

```julia-repl
julia> out = md_isokinetic(sys,Options(iequil=200,nsteps=2_000,ibath=1));
```

The plot of energies as a function of time can be obtained with:
```julia-repl
julia> plot(out,ylim=[-100,100],
           label=["Potential" "Kinetic" "Total" "Temperature"],
           xlabel="step"
       )
```

Note that the total energy is no longer constant during the equilibration period. The potential and kinetic energies should have converged somewhat better than in the simulation without temperature control, although this first simulation is very short.

The temperature can be observed with:
```julia-repl
julia> plot(
           out[:,4],
           label=["Temperature"],
           xlabel="step"
       )
```
Note that it remains practically constant and equal to the target temperature (0.60) during equilibration, but divergences may be observed afterward. If the system is not equilibrated, these divergences can be systematic.

Try different parameters and understand the effect of the equilibration time and the equilibration frequency on the temperature.

A good set of conditions to visualize results is obtained with `ibath=50` and `iequil=5_000`, for `nsteps=20_000`.
```julia-repl
julia> out = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000));
```

Under these conditions, normally, no systematic drift in energies or temperature should be observed after equilibration. Repeat the plots (in the `Julia` prompt, the up arrow accesses previous commands).

The trajectory, `traj.xyz`, can be viewed with `VMD`, as explained earlier.

## 4.2. Complete summary code

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000))
plot(
    out,ylim=[-100,100],
    label=["Potential" "Kinetic" "Total" "Temperature"],
    xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```
