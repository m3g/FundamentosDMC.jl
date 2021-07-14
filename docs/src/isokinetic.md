
# Control de temperatura isocinético

 La función implementada en [md-isokinetic.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-isokinetic.jl) implementa el control de temperatura
isocinético. En este método, las velocidades son escalonadas por un
parámetro $\lambda = \sqrt{T_0/T}$ a intervalos regulares, para
termostatizar el sistema a la temperatura $T_0$. 
Su ejecución demanda la definición de dos parámetros
adicionales: el intervalo de tiempo entre dos escalonamientos de
velocidades, y el tiempo de equilibración. El tiempo de equilibración es
el tiempo dentro del cual los escalonamientos son realizados. El
objetivo debe ser obtener una simulación estable, con energía cinética
media adecuada a la deseada (60 unidades), después de la equilibración.

# 4.1. Control de parámetros y termostatización

El sistema es inicializado de la misma forma que antes, esto es:
```julia-repl
julia> using FundamentosDMC, Plots

julia> sys = System(n=100,sides=[100,100]) 

julia> minimize!(sys);
```

En seguida, vamos a ejecutar la simulación, ahora con termostato isocinético, por `2000` pasos, de los cuales `200` son de equilibración. El termostato es aplicado a cada `ibath` pasos:

```julia-repl
julia> out = md_isokinetic(sys,Options(iequil=200,nsteps=2_000,ibath=1));
```

El gráfico de las energías en función del tiempo puede ser obtenido con:
```julia-repl
julia> plot(out,ylim=[-100,100],
         label=["Potential" "Kinetic" "Total" "Temperature"],
         xlabel="step"
       )
```

Debe notar-se que la energía total no es más constante durante el período de equilibración. Las energía potencial y cinética deben haber convergido un poco mejor que en la simulación sin control de temperatura, aunque esta primera simulación es muy corta. 

La temperatura puede ser observada con:
```julia-repl
julia> plot(
         out[:,4],
         label=["Temperature"],
         xlabel="step"
       )
```
Note que se mantiene prácticamente constante e igual a la temperatura objetivo (0.60) durante la equilibración, pero después divergencias pueden ser observadas. Si el sistema no está equilibrado, estas divergencias pueden ser sistemáticas. 

Pruebe diferentes parámetros, y entienda el efecto del tiempo de equilibración y de la frecuencia de equilibración sobre la temperatura. 

Una buena condición para visualizar los resultados se obtiene con `ibath=50`  y `iequil=5_000`, para `nsteps=20_000`. 
```julia-repl
julia> out = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000));
```

En estas condiciones, normalmente, no se debe observar un desvío sistemático de las energías o de la temperatura después de la equilibración. Repita los gráficos (en el prompt de `Julia`, la flecha para arriba accede a los comandos anteriores).

La trayectoria, `traj.xyz`, puede ser vista con `VMD`, como explicado anteriormente. 

## 4.2. Código completo resumido

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
