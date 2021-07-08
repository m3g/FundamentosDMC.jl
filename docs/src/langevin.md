
# Termostato de Langevin

En la función [md-langevin.jl](https://github.com/m3g/CELFI.jl/blob/master/src/md-langevin.jl) está implementado el control de temperatura
de Langevin.

El termostato de Langevin consiste en la suposición de que dada
partícula real del sistema está inmersa en un fluido de partículas mucho
menores, que desaceleran por fricción y, al mismo tiempo,
ejercen fuerzas aleatorias sobre las partículas reales. La
desaceleración por fricción es proporcional a la velocidad, y los
choques aleatorios causan variaciones instantáneas en las velocidades de
las partículas reales. La trayectoria de una partícula real es, así,
propagada modificando las fuerzas y velocidades. Las fuerzas son
modificadas por la introducción de la fricción,

$$\vec{F}(t) = -\nabla V[\vec{x}(t)] - \lambda \vec{v}(t)$$

Y para que los choques aleatorios de las partículas del baño sean
instantáneos (provoquen cambios instantáneos en las velocidades), una
opción es introducirlos como cambios en la velocidad:

$$v(t+\Delta t) = v(t) + a(t)\Delta t + \sqrt{2\lambda kT \Delta t}\delta(t)$$

Donde $\delta(t)$ es una variable aleatoria Gaussiana con media cero y
desvío estándar 1. La relación entre el coeficiente de fricción, $\lambda$,
que desacelera las partículas, y la intensidad de los choques
estocásticos, $\sqrt{2\lambda kT}$ es resultado del teorema de
fluctuación-disipación, que describe el movimiento Browniano.

# 6.1. Control de parámetros y termalización  

El coeficiente de fricción, $\lambda$, controla el comportamiento de una dinámica de Langevin. En este caso, programa inicializa las velocidades en cero, para que el efecto del termostato sea destacado. Vamos a imprimir la trajectoria con más frequéncia para que se note el efecto de la friccion al início:

Por ejemplo,
```julia-repl
julia> out = md_langevin(sys,Options(lambda=0.1,nsteps=20_000,iprintxyz=5))

```
Ejecute el programa con diversos parámetros, en particular estos:

| $\lambda$ | $\Delta t$ | 
|:---------:|:----------:|  
| 0.001 |  0.05 |
| 0.01  |  0.05 |
| 10.   |  0.05 |
| 10.   |  0.001|
|       |       |

En seguida de cada ejecución, observe los gráficos de energía y las
trayectorias. Discuta si la temperatura llegó al valor deseado (energía
cinética igual a 60 unidades), y si la energía total es en media
constante. Observe el movimiento de las partículas en cada trayectoria. 
La consistencia del termostato depende del paso de integración, en
particular para acoplamientos grandes.  

## 6.2. Código completo resumido

```julia
using CELFI, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_langevin(sys,Options(lambda=0.1,nsteps=20_000,iprintxyz=1))
plot(
  out,ylim=[-100,100],
  label=["Potential" "Kinetic" "Total" "Temperature"],
  xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```
