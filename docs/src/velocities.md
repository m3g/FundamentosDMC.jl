
# Distribución de velocidades

Las velocidades de las partículas en un sistema molecular siguen la distribución de Maxwell-Boltzmann. Vamos a ver si esto efectivamente ocurre en nuestras simulaciones. 

## 10.1. Distribución inicial de velocidades

Un aspecto interesante de ver en las simulaciones es como evoluciona, y como 
converge, la distribución de velocidades de las partículas. En nuestros programas
tenemos trés diferentes maneras de generar la distribución *inicial* de velocidades: con normas distribuidas homogeneamente entre 0 y 1, con una distribución normal de normas, o todas nulas. Esto se puede ajustar con el parámetro `initial_velocities` de `Options`. Por ejemplo, velocidades iniciales con distribuciones distintas pueden ser generadas de la siguiente forma: 

```julia-repl
julia> sys = System(n=100_000,size=[100,100,100])
System{Point2D}:
 Number of particles = 100000
 Box sides = [100.0, 100.0]

julia> v_normal = init_velocities(sys,Options(initial_velocities=:normal))
10000-element Vector{Point2D}:
 [0.04326001565583558, 0.4435409371990999]
 ⋮
 [-0.13167046004485733, -0.8914869714395789]

 julia> v_flat = init_velocities(sys,Options(initial_velocities=:flat))
10000-element Vector{Point2D}:
 [-0.6634658506319203, 0.9721022278467606]
 ⋮
 [1.773444140352424, 0.27389840410745514]

```

Para ver la distribución de velocidades en cada caso, vamos a usar la función `velocity_distribution`, que recibe lee las velocidades grabadas en un archivo de trayectorias (la razón desto va a quedar clara en seguida):

```julia-repl
julia> printxyz(v_normal,sys,"vnormal.xyz")

julia> printxyz(v_flat,sys,"vflat.xyz")

julia> d_normal = velocity_distribution(sys,"vnormal.xyz");
Number of frames read = 1
 Average velocity = [-3.4000000016826124e-8, 4.899999999760096e-8]
 Average velocity norm = 0.9699959035972906
 Average kinetic energy = 0.5998727466282452

julia> d_flat = velocity_distribution(sys,"vflat.xyz");
Number of frames read = 1
 Average velocity = [-4.9000000006849124e-8, -1.2999999991358812e-8]
 Average velocity norm = 0.9506616041200406
 Average kinetic energy = 0.5999232099558015

```

Note que generamos `100_000` puntos, para tener una distribución buen muestreada. 
Las dos distribuciones tienen la misma energia cinética media, claro, que corresponde a la temperatura deseada. Las distribuciones, por otro lado, son muy distintas:

```julia-repl
julia> plot(
         [d_normal,d_flat],
         label=["Normal" "Flat"],
         xlabel="velocity norm",
         ylabel="frequency",
         linewidth=2
       )
```
que resulta en:
```@raw html
<center>
<img src=../figures/velocities.svg>
</center>
```


