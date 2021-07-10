# Comparación entre MC y MD

Vamos a comparar el conjunto de estructuras generados por una dinámica molecular con un conjunto de estructuras generado por una simulación de Monte-Carlo. 

## 8.1. Generando un buen muestreo por MD

Ejecute, nuevamente, el programa `md_langevin` con los
parámetros $\lambda=0.01$, por 100 mil pasos, y salve la trayectoria en un archivo con nombre próprio, como `md.xyz`: 

```julia-repl
julia> md_out = md_langevin(
         sys,
         Options(lambda=0.01,nsteps=100_000,trajectory_file="md.xyz")
       );
```

## 8.2. Generando un buen muestreo de MC

Ejecute una larga simulación de Monte-Carlo (`200_000` pasos)

```julia-repl
julia> mc_out = mc(
         sys,
         Options(alpha=0.05,nsteps=200_000,trajectory_file="mc.xyz")
       );
```

## 8.3. Comparando las energias potenciales

La primera columna de `md_out` tiene la energia potencial en la simulación de dinámica molecular. Vamos a hacer el gráfico de esta energia:

```julia-repl
julia> plot(md_out[:,1],label="Potential energy- MD",xlabel="step")
```

Y vamos a añadir al mismo gráfico la energia potencial obtenida en la simulación de Monte-Carlo, que está en `mc_out` (note el `!` al final de `plot!`, indicando que el gráfico anterior va a ser modificado):

```julia-repl
julia> plot!(mc_out,label="Potential energy- MC",xlabel="step")
```

Note la similaridad, o diferencia entre los dos gráficos. Acuerde-se que en ninguna de estas simulaciones controlamos explicitamente la energia potencial. 

## 8.4. Función de distribución radial 

Vamos a comparar la estructura media obtenida usando MD con la
estructura media obtenida con MC. Para eso vamos a usar la función
$g(r)$, que se llama ``función de distribución radial''. 
Esta función mide la probabilidad de encontrar una partícula a
una distancia $r$ de otra partícula en el sistema real, en relación a
esa misma probabilidad si no hubiese ninguna interacción entre las
partículas. 

En nuestro caso bi-dimensional, el número de partículas por unidad de
volumen es $\rho=n/A$, donde $n=100$ es el número de partículas y $A=100^2$
es el área total del sistema simulado. El número de partículas esperado
en un intervalo de distancias entre $r$ y $r+\Delta r$ de cada partícula
es, por lo tanto $n(r)=\rho A(r)$, donde $A(r)$ es el área de una cáscara
circular de radio menor $r$ y radio mayor $r+\Delta r$:
```@raw html
<center>
<img src=../figures/area.svg>
</center>
```
Vemos que $A(r)=\pi (r+\Delta r)^2 - \pi r^2 \approx 2\pi r\Delta r$.
De este modo, el número de partículas esperado, en media, seria de 
$n(r)=2\pi r\Delta r\rho$, si no hubiese interacciones. 

Las interacciones hacen con que el número de partículas en cada
distancia sea diferente de una distribución homogénea. Si hay
interacciones favorables, por ejemplo, la probabilidad de encontrar dos
partículas próximas es mayor. Esta distribución de partículas es uno de
los parámetros estructurales más importantes.

## 8.5. Cálculo de $g(r)$

El programa [radial_distribution.jl](https://github.com/m3g/FUndamentosDMC.jl/blob/master/src/radial_distribution.jl)   calcula, a partir de una trayectoria, la función
$g(r)=n'(r)/n(r)$, done $n(r)$ esta definido anteriormente, y $n'(r)$ es
el número medio de partículas efectivamente observado entre $r$ y $r+\Delta r$
en la simulación. 

Para calcular la función de distribución radial de la simulación de dinámica molecular, ejecute:
```julia-repl    
julia> rmd, gmd = radial_distribution(sys,"md.xyz")
```

`rmc` y `rmd` son las distancias y las densidades relativas de partículas a cada distancia, en relación a la densidad del sistema. 

En seguida, obtengamos la $g(r)$ de la simulación de Monte-Carlo: 
```julia-repl    
julia> rmc, gmc = radial_distribution(sys,"mc.xyz")
```

Cree en gráfico de la función de distribución de la simulación de dinámica molecular, con:
```julia-repl
julia> plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
```

Note el aumento de densidad local en distancias cortas, que resulta de las interacciones favorables entre las partículas. Note también que la densidad relativa tiende a `1.0` en distancias grandes, cuando la correlación de la posición entre partículas no és mas importante.  

Agregue al mismo gráfico la función de distribución obtenida a partir de la simulación de Monte-Carlo:  
```julia-repl
julia> plot!(rmc,gmc,xlabel="r",ylabel="g(r)",label="MC")
```
 Compare. Las simulaciones, con sus naturalezas
totalmente distintas, muestrearon las mismas estructuras?  

## 8.2. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)

md_out = md_langevin(
  sys,
  Options(lambda=0.01,nsteps=100_000,trajectory_file="md.xyz")
)
mc_out = mc(
  sys,
  Options(alpha=0.05,nsteps=200_000,trajectory_file="mc.xyz")
)
plot(md_out[:,1],label="Potential energy- MD",xlabel="step")
plot!(mc_out,label="Potential energy- MC",xlabel="step")
savefig("potential_energy.pdf")

rmd, gmd = radial_distribution(sys,"md.xyz")
rmc, gmc = radial_distribution(sys,"mc.xyz")
plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
plot!(rmc,gmc,xlabel="r",ylabel="g(r)",label="MC")
savefig("gr.pdf")
```

