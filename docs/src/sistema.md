
# Sistema simulado

La simulación es de un fluido de 100 partículas (mono-atómicas) que
interactúan por un
potencial de Lennard-Jones, en un sistema bi-dimensional, periódico.

$V = 4\epsilon \left( \frac{\sigma^{12}}{r^{12}} - \frac{\sigma^6}{r^6} \right)$

Abra el archivo [potential.jl](https://github.com/m3g/CELFI.jl/blob/master/src/potential.jl) y entienda la implementación del
cálculo de la energía potencial. Note que el cálculo depende de 3
parámetros: $\epsilon$, $\sigma$, y el tamaño del sistema periódico. Los
parámetros están definidos en la estructura de datos `opt`, de entrada (veremos más tarde como usarla). 

El archivo [forces.jl](https://github.com/m3g/CELFI.jl/blob/master/src/forces.jl) contiene el cálculo de las fuerzas (el gradiente
del potencial), y el archivo [kinetic.jl](https://github.com/m3g/CELFI.jl/blob/master/src/kinetic.jl) contiene el cálculo
de la energía cinética. Como el sistema usa condiciones periódicas de
contorno, las coordenadas tienen que siempre ser calculadas en relación
a la imagen mínima. El cálculo de la imagen mínima está implementado en
el archivo [image.jl](https://github.com/m3g/CELFI.jl/blob/master/src/image.jl). Es interesante entender la
implementación de cada una de estas funciones, que son comunes a todos
los métodos que vamos a describir. 

## 2.1. Parámetros y opciones de la simulación

Para inciar los trabajos, abra una sección de `Julia`, y dé el comando:
```julia-repl
julia> using CELFI
```

Los parámetros de las simulaciones son controlados en la inicialización de la estructure `Options`, por ejemplo:

```julia-repl
julia> opt = Options(sides=[100,100],dt=0.01)
Options{Point2D}
-------------------
Simulation options:
-------------------
dt = 0.01
nsteps = 2000
sides = [100.0, 100.0]
eps = 1.0
sig = 2.0
kavg_target = 0.6
ibath = 1
printxyz = true
printvel = false
iprint = 1
iprintxyz = 2
trajectory_file = traj.xyz
energies_file = energies.dat
velocities_file = velocities.dat

julia>
```

En este caso, ajustamos en tamaño del sistema y el paso de tiempo manualmente, y mantuvimos todas las otras opciones con valores default. Cada uno de estos paráemetros será discutido oportunamente. Note que definen el tamaño, campo de fuerza ($\epsilon$ y $\sigma$), energía cinética (temperatura), y los nombres de los archivos de salida. 

## 2.2. Coordenadas iniciales

La coordenadas iniciales pueden ser creadas aleatoriamente, usando: 
```julia-repl
julia> x = [ opt.sides .* rand(Point2D) for i in 1:100 ]
100-element Vector{Point}:
 [18.36579648764145, 7.711401822973363]
 [41.784092301135665, 45.61852102711508]
 [23.850299728474454, 63.797752122286425]
 ⋮
 [92.5679156243071, 39.272476774702206]
 [26.845528447086274, 92.88216539818639]

```
que genera `100` puntos en 2 dimensiones, aleatórios, con coordenadas entre `[0,0]` y `opt.sides = [100,100]`, en este caso. `Point2D` es un tipo de variable que representa un punto en dos dimensiones. Mas adelante vamos a ver que todo el código es genérico, y podemos hacer simulaciones en 3 dimensiones apenas modificando el tipo de variable asociado. 

!!! note
    El punto `.` en `.*` indica que es una multiplicación componente-a-componente, de cada componente del vector `opt.sides` por cada componente del vector associado a cada punto. 

## 2.3. Minimización de la energia

En seguida, el punto inicial va a ser modificado usando 
el Método del Gradiente para minimizar la energía. El
método consiste en mover las partículas según la aproximación de Taylor
de orden uno, en la dirección de descenso de energía:

$\vec{x}_{i+1} = \vec{x}_i - \nabla V(\vec{x}_i) \Delta x$

Si la energía en el punto $\vec{x}_{i+1}$ es menor que la energía en el
punto $\vec{x}_i$, se acepta el punto $\vec{x}_{i+1}$ y el proceso es
repetido. Si no, $\Delta x$ es disminuido ($\Delta x = \Delta x / 2$), y
un nuevo punto $\vec{x}_{i+1}$ es calculado. Como la aproximación debe
ser una buena aproximación en las cercanias del punto corriente ($\vec{x}_i$), un
gradiente negativo garante que la función disminuye para $\Delta x$
suficientemente pequeño. El proceso es interrumpido cuando la norma del
gradiente es pequeña, o cuando demasiados puntos fueron testados. En mecánica, $-\nabla V = \vec{F}$, entonces la función que calcula el gradiente es la misma que calcula las fuerzas en la simulación. Abra el archivo [minimize.jl](https://github.com/m3g/CELFI.jl/blob/master/src/minimize.jl) para
discutir como se crea el punto inicial. 

Antes de ejecutar la minimización de energia, vamos a copiar el punto inicial, para comparación:
```julia-repl
julia> x0 = copy(x)
```

En seguida, minimizamos la energia con la función `minimize!`:
```julia-repl
julia> minimize!(x,opt)
```

En `Julia` es una convención que las funciones que modifican sus argumentos terminan en `!`. En este caso, la función va a modificar las posiciones, `x`, de las partículas. 

## 2.4. Temperatura

La temperatura del sistema es un parámetro también definido internamente
en el programa (puede ser modificado a gusto, pero no lo haremos). La
temperatura se define a partir energía cinética media asociada a cada
grado de libertad de movimiento del sistema. En el caso que todos los
movimientos pueden ser escritos como translaciones, la definición es

$\frac{1}{2}kT = \left< \frac{1}{2} m v_x^2\right>$

donde la media, hecha sobre $v_x$ aqui, es equivalente si hecha sobre
cualquier otro grado de libertad de translación. En un sistema
tridimensional, por lo tanto, 

$\left<\frac{1}{2}m |\vec{v}|^2 \right> = 
\left<\frac{1}{2}m \left(v_x^2 + v_y^2 + v_z^2\right) \right> = 
3\left< \frac{1}{2} m v_x^2 \right> = \frac{3}{2}kT$

que es el resultado usual.

Nuestras simulaciones son de un sistema bi-dimensional. En este caso,

$\left< \frac{1}{2}m |\vec{v}|^2 \right> = 
\left< \frac{1}{2}m \left(v_x^2 + v_y^2\right)\right> =
2\left< \frac{1}{2}m v_x^2 \right> = kT$

En los códigos de dinámica molecular, la definición de temperatura se
da, así, por la definición de la energía cinética media o, en este caso,
por $kT$. En el código de Monte-Carlo la definición de temperatura se da
por la tasa de aceptación, con la misma definición. 

En todos los códigos fue escogido que se objetiva simular el sistema a
la temperatura que corresponde a $kT = 0.6$ unidades. Los sistemas
simulados tiene 100 partículas, por lo tanto la energía cinética media
es $100kT=60$ unidades. Las velocidades iniciales van a ser generadas aleatoriamente al princípio de la simulación. 






