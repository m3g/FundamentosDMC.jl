
# Sistema simulado

La simulación es de un fluido de 100 partículas (mono-atómicas) que
interactúan por un
potencial de Lennard-Jones, en un sistema bi-dimensional, periódico.

$
V = 4\epsilon \left( \frac{\sigma^{12}}{r^{12}} - \frac{\sigma^6}{r^6} \right)
$

[teste](t)

Abra el archivo [potential.jl](https://github.com/m3g/CELFI.jl/blob/master/src/potential.jl) y entienda la implementación del
cálculo de la energía potencial. Note que el cálculo depende de 3
parámetros: $\epsilon$, $\sigma$, y el tamaño del sistema periódico. Los
parámetros están definidos en la estructura de datos `opt`, de entrada (veremos más tarde como usarla). 

El archivo [forces.jl](https://github.com/m3g/CELFI.jl/blob/master/src/forces.jl) contiene el cálculo de las fuerzas (el gradiente
del potencial), y el archivo [kinetic.jl](https://github.com/m3g/CELFI.jl/blob/master/src/kinetic.jl) contiene el cálculo
de la energía cinética. Como el sistema usa condiciones periódicas de
contorno, las coordenadas tienen que siempre ser calculadas en relación
a la imagen mínima. El cálculo de la imagen mínima está implementado en
el archivo [image.jl](https://github.com/m3g/CELFI.jl/blob/master/src/kinetic.jl). Es interesante entender la
implementación de cada una de estas funciones, que son comunes a todos
los métodos que vamos a describir. 

## Temperatura

La temperatura del sistema es un parámetro también definido internamente
en el programa (puede ser modificado a gusto, pero no lo haremos). La
temperatura se define a partir energía cinética media asociada a cada
grado de libertad de movimiento del sistema. En el caso que todos los
movimientos pueden ser escritos como translaciones, la definición es

$\frac{1}{2}kT = \left< \frac{1}{2} m v_x^2\right>$

donde la media, hecha sobre $v_x$ aqui, es equivalente si hecha sobre
cualquier otro grado de libertad de translación. En un sistema
tridimensional, por lo tanto, 

$
\left<\frac{1}{2}m |\vec{v}|^2 \right> = 
\left<\frac{1}{2}m \left(v_x^2 + v_y^2 + v_z^2\right) \right> = 
3\left< \frac{1}{2} m v_x^2 \right> = \frac{3}{2}kT
$

que es el resultado usual.

Nuestras simulaciones son de un sistema bi-dimensional. En este caso,

$
\left< \frac{1}{2}m |\vec{v}|^2 \right> = 
\left< \frac{1}{2}m \left(v_x^2 + v_y^2\right)\right> =
2\left< \frac{1}{2}m v_x^2 \right> = kT 
$

En los códigos de dinámica molecular, la definición de temperatura se
da, así, por la definición de la energía cinética media o, en este caso,
por $kT$. En el código de Monte-Carlo la definición de temperatura se da
por la tasa de aceptación, con la misma definición. 

En todos los códigos fue escogido que se objetiva simular el sistema a
la temperatura que corresponde a $kT = 0.6$ unidades. Los sistemas
simulados tiene 100 partículas, por lo tanto la energía cinética media
es $100kT=60$ unidades.



