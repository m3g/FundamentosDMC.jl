# Simulación de Dinámica Molecular Microcanónica

Abra el archivo [md-simple.jl](https://github.com/m3g/CELFI.jl/blob/master/src/md-simple.jl), que contiene el algoritmo de
simulación. La simulación empieza con velocidades aleatorias, ajustadas
para la media termodinámica deseada de 0.6 unidades/átomo. A esta
energía cinética media le llamaremos ``temperatura''.
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

La simulación no tiene control de temperatura o de presión. Es una propagación de la trajectória según las leyes de Newton, que deberían conservar la energia. A esto se le llama una simulación "microcanónica", o "NVE" (que conserva, en princípio, el número de partículas, el volumen y la energia total).

## 3.1. Paso de integración

Para realizar una MD simple, con un paso de integración de `dt=1.0`, ejecute le comando:
```julia-repl
julia> md(x,Options(dt=1.0))

```
En princípio, está previso realizar 2000 pasos de integración de las equaciones 
de movimimiento.  Pruebe pasos de integración entre `1.0` y `0.01`. 
Note que pasa con la energía. Note que pasa con la energía
cinética media, la cual fue inicializada en 0.6 unidades/átomo. Discuta
la elección del paso de integración, y los valores de energía cinética
obtenidos. Las simulaciones que siguen van a usar un paso de integración
`dt = 0.05`.

Al conseguir una simulación estable hasta
el fin, observe el gráfico de energía en función del tiempo, usando el
comando:
\command{xmgrace -nxy energies.dat}
Observe y trate de entender las amplitudes de las oscilaciones de las
energías cinética y potencial, y las amplitudes de las
oscilaciones de la energía total. A que se deben cada una de las
oscilaciones? Observe como estas oscilaciones dependen del paso de
integración.

Por fin, visualice la trayectoria, usando
\command{vmd -e view.vmd}
Dentro de VMD, dé el comando 
\command{pbc set \{ 100. 100. 100. \} -all}
y represente explícitamente el sistema periódico eligiendo {\tt
+X;+Y;-X;-Y} en
\command{Graphics$\to$Representations$\to$Periodic}


