
# Sistema simulado

La simulación es de un fluido de 100 partículas (mono-atómicas) que
interactúan por un
potencial de Lennard-Jones, en un sistema bi-dimensional, periódico.

$V = 4\varepsilon \left( \frac{\sigma^{12}}{r^{12}} - \frac{\sigma^6}{r^6} \right)$

Abra el archivo [potential.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/potential.jl) y entienda la implementación del
cálculo de la energía potencial. Note que el cálculo depende de 3
parámetros: $\varepsilon$, $\sigma$, y el tamaño del sistema periódico. Los
parámetros están definidos en la estructura de datos `opt`, de entrada (veremos más tarde como usarla). 

El archivo [forces.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/forces.jl) contiene el cálculo de las fuerzas (el gradiente
del potencial), y el archivo [kinetic.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/kinetic.jl) contiene el cálculo
de la energía cinética. Como el sistema usa condiciones periódicas de
contorno, las coordenadas tienen que siempre ser calculadas en relación
a la imagen mínima. El cálculo de la imagen mínima está implementado en
el archivo [image.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/image.jl). Es interesante entender la
implementación de cada una de estas funciones, que son comunes a todos
los métodos que vamos a describir. 

## 2.1. Coordenadas iniciales

Para iniciar los trabajos, abra una sección de `Julia`, y dé el comando:
```julia-repl
julia> using FundamentosDMC
```
El sistema inicial puede ser creado aleatoriamente, usando: 
```julia-repl
julia> sys = System(n=100,sides=[100,100]) 
System{Point2D}:
 Number of particles = 100
 Box sides = [100.0, 100.0]

```
que genera `100` puntos en 2 dimensiones, aleatorios, con coordenadas entre `[-50,-50]` y `[50,50]`, en este caso. `Point2D` es un tipo de variable que representa un punto en dos dimensiones. Mas adelante vamos a ver que todo el código es genérico, y podemos hacer simulaciones en 3 dimensiones apenas modificando el tipo de variable asociado. 

## 2.1. Parámetros y opciones de la simulación

Los parámetros de las simulaciones son controlados en la estructura `Options`, por ejemplo, para ajustar el paso de tiempo, pasamos el parámetro `dt` a la estructura. Esto puede ser hecho en la llamada de las funciones de simulación, como veremos.

```julia-repl
julia> Options(dt=0.1)
Options
  dt: Float64 0.1
  nsteps: Int64 2000
  eps: Float64 1.0
  sig: Float64 2.0
  initial_velocities: Symbol normal
  kT: Float64 0.6
  ibath: Int64 10
  iequil: Int64 10
  tau: Int64 10
  lambda: Float64 0.1
  alpha: Float64 0.1
  printxyz: Bool true
  iprintxyz: Int64 10
  iprint: Int64 20
  minprint: Bool false
  trajectory_file: String "traj.xyz"
  velocity_file: String "vel.xyz"
```

En este caso, ajustamos el paso de tiempo manualmente, y mantuvimos todas las otras opciones con valores default. Cada uno de estos parámetros será discutido oportunamente. Note que definen el tamaño, campo de fuerza ($\varepsilon$ y $\sigma$), $kT$ (temperatura), y los nombres de los archivos de salida. 

## 2.3. Minimización de la energía

En seguida, el punto inicial va a ser modificado usando 
el Método del Gradiente para minimizar la energía. El
método consiste en mover las partículas según la aproximación de Taylor
de orden uno, en la dirección de descenso de energía:

$\vec{x}_{i+1} = \vec{x}_i - \nabla V(\vec{x}_i) \Delta x$

Si la energía en el punto $\vec{x}_{i+1}$ es menor que la energia en el
punto $\vec{x}_i$, se acepta el punto $\vec{x}_{i+1}$ y el proceso es
repetido. Si no, $\Delta x$ es disminuido ($\Delta x = \Delta x / 2$), y
un nuevo punto $\vec{x}_{i+1}$ es calculado. Como la aproximación debe
ser una buena aproximación en las cercanias del punto corriente ($\vec{x}_i$), un
gradiente negativo garante que la función disminuye para $\Delta x$
suficientemente pequeño. El proceso es interrumpido cuando la norma del
gradiente es pequeña, o cuando demasiados puntos fueron testados. En mecánica, $-\nabla V = \vec{F}$, entonces la función que calcula el gradiente es la misma que calcula las fuerzas en la simulación. Abra el archivo [minimize.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/minimize.jl) para
discutir como se crea el punto inicial. 

Antes de ejecutar la minimización de energía, vamos a copiar el punto inicial, para comparación:
```julia-repl
julia> x0 = copy(sys.x0)
```

En seguida, minimizamos la energía con la función `minimize!`:
```julia-repl
julia> minimize!(sys)
Energy before minimization: 38322.72337856496
Energy after minimization: -74.15646912098042
```

!!! note
    En `Julia` es una convención que las funciones que modifican sus argumentos terminan en `!`. En este caso, la función va a modificar las posiciones, `sys.x0`, de las partículas. 

Si desea ver la progresión de las energías, active la opción `minprint`, con `minimize!(sys,Options(minprint=true))`.

Podemos ver rápidamente que ocurrió con las partículas, colocando-las en un gráfico. Primero, generamos un gráfico de los puntos antes de la minimización:
```julia-repl
julia> using Plots

julia> scatter(Tuple.(x0))
```
Los puntos deben estar aleatoriamente distribuidos, y en particular algunos puntos deben estar muy cercanos a los otros, lo que genera potenciales muy repulsivos.

En seguida, hacemos el gráfico del punto con energía mínima obtenido:
```julia-repl
julia> scatter(Tuple.(sys.x0))
```
y notará que los puntos ahora tienen una nueva disposición: hay puntos formando clusteres, porque el potencial de Lennard-Jones es atractivo en distancias largas. Pero no hay más puntos muy cercanos generando repulsiones muy grandes.

Este punto inicial de energía mínima será usado en nuestras simulaciones. 

!!! note
    Los gráficos pueden ser salvos como figuras con `savefig("plot.pdf")`, por ejemplo.

## 2.4. Temperatura

La temperatura del sistema es un parámetro también definido internamente
en el programa (puede ser modificado a gusto, pero no lo haremos). La
temperatura se define a partir de la energía cinética media asociada a cada
grado de libertad de movimiento del sistema. En el caso que todos los
movimientos pueden ser escritos como traslaciones, la definición es

$\frac{1}{2}kT = \left< \frac{1}{2} m v_x^2\right>$

donde la media, hecha sobre $v_x$ aqui, es equivalente si hecha sobre
cualquier otro grado de libertad de traslación. En un sistema
tridimensional, por lo tanto, 

$\left<\frac{1}{2}m |\vec{v}|^2 \right> = 
\left<\frac{1}{2}m \left(v_x^2 + v_y^2 + v_z^2\right) \right> = 
3\left< \frac{1}{2} m v_x^2 \right> = \frac{3}{2}kT$

que es el resultado usual.

Nuestras simulaciones son de un sistema bi-dimensional. En este caso,

$\left< \frac{1}{2}m |\vec{v}|^2 \right> = 
\left< \frac{1}{2}m \left(v_x^2 + v_y^2\right)\right> =
2\left< \frac{1}{2}m v_x^2 \right> = kT$

En todos los códigos fue escogido que se objetiva simular el sistema a
la temperatura que corresponde a $kT = 0.6$ unidades. Los sistemas
simulados tienen 100 partículas y son bidimensionales, por lo tanto la energía cinética media es de $100kT=60$ unidades. Las velocidades iniciales van a ser generadas aleatoriamente al principio de la simulación. 

## 2.5. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
x0 = copy(sys.x0)
minimize!(sys)
scatter(Tuple.(x0))
scatter(Tuple.(sys.x0))
```

Observación: el comando `Tuple.(x)` convierte el vector de vectores en un vector de pares (tuplas), que es correctamente interpretado por `Plots` como una única serie de puntos.



