
# Control de temperatura isocinético

 La función implementada en [md-isokinetic.jl](https://github.com/m3g/CELFI.jl/blob/master/src/md-simple.jl) implementa el control de temperatura
isocinético. En este método, las velocidades son escalonadas por un
parámetro $\lambda = \sqrt{T_0/T}$ a intervalos regulares, para
termostatizar el sistema a la temperatura $T_0$. 
Su ejecución demanda la definición de dos parámetros
adicionales: el intervalo de tiempo entre dos escalonamientos de
velocidades, y el tiempo de equilibración. El tiempo de equilibración es
el tiempo dentro del cual los escalonamientos son realizados. El
objetivo debe ser obtener una simulación estable, con energía cinética
media adecuada a la deseada (60 unidades), después de la equilibración.

El sistema es inicializado de la misma forma que antes, esto es:
```julia-repl
julia> using CELFI, Plots

julia> opt = Options(side=[100,100]);

julia> x = [ opt.sides .* rand(Point2D) for i in 1:100 ];

julia> minimize!(x,opt);
```

En seguida, vamos a ejecutar la simulación, ahora con termostato isocinético, por `2000` pasos, de los cuales `200` son de equilibración. El termostato es aplicado a cada `ibath` pasos:

```julia
julia> log = md_isokinetic(x,Options(iequil=200,nsteps=2_000,ibath=20))
```

El gráfico de las energias en función del tiempo puede ser obtenido con:
```julia-repl
julia> plot(log,ylim=[-100,100],
         label=["Potential" "Kinetic" "Total" "Temperature"],
         xlabel="step"
       )

```

Debe notar-se que la energia total no es más constante dureante el periódo de equilibración. Las energia potencial y cinética deben haber convergido un poco mejor que en la simulación sin control de temperatura, aunque esta primera simulación es muy corta. 

La temperatura puede ser observada con:
```julia-repl
julia> plot(
         log[:,4],
         label=["Temperature"],
         xlabel="step"
       )
```
Note que se mantiene practicamente constante e igual a la temperatura objetivo (0.60) durante la equilibración, pero después divergéncias pueden ser observadas. Si el sistema no está equilibrado, estas divergéncias pueden ser sistemáticas. 

Pruebe diferentes parámetros, y entienda el efecto del tiempo de equilibración y de la frequéncia de equilibración sobre la temperatura. 

Una buena condición para visualizar los resultados se obtiene con `ibath=20`  y
`iequil=20_000`, para `nsteps=20_000`. 
```julia
julia> log = md_isokinetic(x,Options(iequil=2000,ibath=20,nsteps=20_000))
```

En estas condiciones, normalmente, no se debe observar un desvio sistemático de las energias o de la temperatura después de la equilibración. 

La trajectoria, `traj.xyz`, puede ser vista con `VMD`, como explicado anteriormente. 

## Código completo resumido

```julia
using CELFI, Plots
opt = Options(sides=[100,100])
x = [ opt.sides .* rand(Point2D) for i in 1:100 ]
minimize!(x,opt)
log = md_isokinetic(x,Options(iequil=2000,ibath=20,nsteps=20_000))

plot(
  log,ylim=[-100,100],
  label=["Potential" "Kinetic" "Total" "Temperature"],
  xlabel="step"
)

plot(log[:,4],label="Temperature",xlabel="step")
```
