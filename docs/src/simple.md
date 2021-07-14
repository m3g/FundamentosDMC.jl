# Simulación de Dinámica Molecular Microcanónica

Abra el archivo [md-simple.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-simple.jl), que contiene el algoritmo de
simulación. La simulación empieza con velocidades aleatorias, ajustadas
para la media termodinámica deseada de 0.6 unidades/átomo ($kT=0.6$ en un sistema bidimensional). A esta
energía cinética media le llamaremos "temperatura".
El algoritmo de integración es Velocity-Verlet, que consiste en propagar
la posiciones con

$\vec{x}(t+\Delta t) = \vec{x}(t) + \vec{v}(t)\Delta t + \frac{1}{2}\vec{a}(t)\Delta t^2$

siendo $\vec{a}(t)=\vec{F}(t)/m$, donde $\vec{F}(t)$ es la fuerza en el tiempo corriente. 
La fuerza en seguida es calculada en un tiempo posterior de tiempo con

$\vec{F}(t+\Delta t) = -\nabla V\left[\vec{x}(t)\right]$

y entonces las velocidades en el instante siguiente son calculadas con

$\vec{v}(t+\Delta t) = \vec{v}(t) +
\frac{1}{2}\left[
\frac{\vec{F}(t)}{m}+\frac{\vec{F}(t+\Delta t)}{m}\right]$

completando el ciclo. En este ejemplo las masas son consideradas
unitarias, para simplificar. La simulación es ejecutada por `nsteps`
pasos, con paso de integración $\Delta t$, este siendo un parámetro de
entrada, `dt`, definido en `Options`. 

La simulación no tiene control de temperatura o de presión. Es una propagación de la trayectoria según las leyes de Newton, que deberían conservar la energía. A esto se le llama una simulación "microcanónica", o "NVE" (que conserva, en principio, el número de partículas, el volumen y la energía total).

## 3.1. Paso de integración

Para realizar una MD simple, con un paso de integración de `dt=1.0`, ejecute el comando:
```julia-repl
julia> out = md(sys,Options(dt=0.1));

```
En principio, está previsto realizar 2000 pasos de integración de las ecuaciones 
de movimiento.  Pruebe pasos de integración entre `1.0` y `0.01`. 
Note que pasa con la energía. Note que pasa con la energía
cinética media, la cual fue inicializada en 0.6 unidades/átomo. Discuta
la elección del paso de integración, y los valores de energía cinética
obtenidos. Las simulaciones que siguen van a usar un paso de integración
`dt = 0.05`.

Es posible controlar la frecuencia de impresión y el número de puntos salvos en el archivo de trayectoria con las opciones `iprint` y `iprintxyz`:
```julia-repl
julia> out = md(sys,Options(dt=0.1,iprint=1,iprintxyz=5))
```
El número total de pasos se controla con el parámetro `nsteps`.

La variable `out` que sale de la simulación es una matriz con las energías y la temperatura en cada paso de la simulación. Es probable que la simulación "explote" con pasos de tiempo grandes. Para visualizar este proceso, podemos hacer:
```julia-repl
julia> using Plots

julia> plot(
         out,ylim=[-100,100],
         label=["Potential" "Kinetic" "Total" "Temperature"],
         xlabel="step"
       )
```

Y obtendremos un gráfico similar a:
```@raw html
<center> 
<img src="../figures/exploded.png">
</center>
```

Para pasos de tiempo menores la simulación debe conseguir llegar hasta el fin. Podemos ver el resultado nuevamente, y debe ser algo similar a: 
```@raw html
<center> 
<img src="../figures/not_exploded.png">
</center>
```

Observe y trate de entender las amplitudes de las oscilaciones de las
energías cinética y potencial, y las amplitudes de las
oscilaciones de la energía total. A que se deben cada una de las
oscilaciones? Observe como estas oscilaciones dependen del paso de
integración.

## 3.2. Visualización de la trayectoria

Por fin, abra la trayectoria usando VMD. No es necesario salir de la sección de `Julia`. Al apretar `;` (punto y coma) aparecerá un prompt `shell>`, desde el cual 
es posible ejecutar a `VMD`, si este está instalado correctamente y disponible en el `path`:

```julia-repl
shell> vmd traj.xyz
```
Dentro de VMD, elija la representación de `VDW`, en 
```
Graphics -> Representations -> Drawing Method -> VDW
```
y dé el comando 
```
vmd> pbc set { 100. 100. 100. } -all
``` 
para indicar la periodicidad del sistema. 
Para representar explícitamente el sistema periódico, elija `+X;+Y;-X;-Y` en
```
Graphics -> Representations -> Periodic
```

Para salir de `VMD` use el comando `exit`, y para volver al prompt de `Julia` desde `shell>`, use `backspace`. 

## 3.3. Código completo resumido

```julia
using FundamentosDMC, Plots

sys = System(n=100)
minimize!(sys)
out = md(sys,Options(dt=0.05))

plot(out,ylim=[-100,100],
  label=["Potential" "Kinetic" "Total" "Temperature"],
  xlabel="step"
)
```

Es generado el archivo `traj.xyz` que puede ser visualizado en `VMD`. 




