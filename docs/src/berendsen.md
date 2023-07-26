
# Control de temperatura de Berendsen

En la función [md-berendsen.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-berendsen.jl) está implementado el control de temperatura
de Berendsen. Este método también es un método basado en el
escalonamiento de velocidades, pero es más suave. Las velocidades son
escalonadas por

$$\lambda = \left[  
1 + \frac{\Delta t}{\tau} \left(
\frac{T_0}{T(t)} -1
\right)
\right]^{1/2}$$

en donde $\Delta t$ es el paso de integración y $\tau$ es un parámetro que
define la velocidad con que el escalonamiento es realizado. El
escalonamiento es más suave y más lento. 

# 5.1. Control de parámetros y termostatización 

El parámetro $\tau$ se ajusta con la opción `tau` de `Options`. Por ejemplo:

```julia-repl
julia> sys = System(n=100,sides=[100,100])

julia> minimize!(sys)

julia> out = md_berendsen(sys,Options(tau=50,iequil=500,nsteps=20_000));
```

Pruebe diferentes parámetros, con `20_000` pasos de simulación. Entre los cuales, estos: 

| $\tau$ |  $i_{\mathrm{equil}}$ | 
|:------:|:---------------------:|
|   50   |  500 |
|   50   |  1500|
|  300   |  1500|
|  300   |  3000|
|        |      |

Observe los gráficos de energía resultantes, usando los mismos comandos de antes:
```julia-repl
julia> using Plots

julia> plot(
           out,ylim=[-100,100],
           label=["Potential" "Kinetic" "Total" "Temperature"],
           xlabel="step"
       )
```

Observe la suavidad, o no, de la curva de energía total. Observe si la
energía cinética se aproximó de la energía media deseada ($kT=60$).

## 5.2. Código completo resumido

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

