# Ejecución en 3D

Todos los códigos aqui pueden generar simulaciones tri-dimensionales. En realidad los únicos códigos que tienen alguna 
dimensionalidad definida son los que escriben el archivo de salida (`printxyz`) y el que calcula la función de distribución
radial (`radial_distribution`). Todos los otros códigos son genéricos, y lo que hace depende apenas del tipo de variable de entrada, o sea, de la dimensión de la caja y de lo puntos. A lo largo de los códigos.  

## 9.1. Inicializando el sistema con 3 dimensiones:

Basta crear el sistema en un caja tri-dimensional:

```julia-repl
julia> sys = System(n=100,sides=[50,50,50])
System{Point3D}:
 Number of particles = 100
 Box sides = [50.0, 50.0, 50.0]
```
Note que ahora los puntos son del tipo `Point3D`. Efectivamente, las coordenadas, `sys.x0` son tri-dimensionales:

```julia-repl
julia> sys.x0
100-element Vector{Point3D}:
 [1.7209592667160578, 2.6345231665739135, 1.9123751433362468]
 [0.8615941909462954, -5.3152138387533565, -1.6775684653805492]
 [5.492878544337906, -7.041239125665051, -4.395671467880198]
⋮
 [-2.229487380572972, -2.4612177799986634, 8.742176399685668]
 [3.9099597073208123, 7.797739759305074, 4.3056981321251655]
```

La energia puede entonces ser mininimizada como antes:
```julia-repl
julia> minimize!(sys)
Energy before minimization: 1.8502746330219281e6
Energy after minimization: -451.7109718153855
```

## 9.2. Ejecutando y visualizando la simulación

La simulación, cualquiera que sea, puede ser ejecutada normalmente. Los códigos ejecutados son exactamente los mismos. Por ejemplo, una simulación con dinámica de Langevin puede ser ejecutada con:

```julia-repl
julia> md_out = md_langevin(
         sys,
         Options(lambda=0.01,nsteps=20_000,trajectory_file="md.xyz")
       )
```

La simulación pude ser visualizada en `VMD`, con:
```
shell> vmd md.xyz
```

Y podemos visualizar la función de distribucion radial, con

```julia-repl    
julia> rmd, gmd = radial_distribution(sys,"md.xyz")

julia> plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")

```

## 9.3. Código completo resumido

```julia
using FundamentosDMC, Plots

sys = System(n=100,sides=[50,50,50])
minimize!(sys)

md_out = md_langevin(
  sys,
  Options(lambda=0.01,nsteps=100_000,trajectory_file="md.xyz")
)
rmd, gmd = radial_distribution(sys,"md.xyz")
plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
```

