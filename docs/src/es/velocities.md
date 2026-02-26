
# Distribución de velocidades

Las velocidades de las partículas en un sistema molecular (tridimensional) siguen la distribución de Maxwell-Boltzmann:

$$f(v) ~dv= \left(\frac{m}{2 \pi kT}\right)^{3/2}\, 4\pi v^2 e^{ -\frac{mv^2}{2kT}} ~ dv$$

Esencialmente, esta ecuación dice que la probabilidad de que una partícula tenga una velocidad en el intervalo $v\pm dv$ es proporcional a $v^2 e^{-v^2/kT}$. Esto se debe a dos motivos: primero, la probabilidad de que una partícula tenga una energía $E$ es dada por la distribución de Boltzmann, $e^{-E/kT}$. Como la energía cinética és $mv^2/2$, el cuadrado de la velocidad aparece en el exponente. Segundo, la *degeneración* aumenta con $v^2$ con el aumento de la velocidad, porque el volumen del espacio que contiene los vectores de longitud $v$ aumenta con $v^2$ en tres dimensiones. Esto genera el término pre-exponencial proporcional a $v^2$. 

Vamos a ver si esto efectivamente ocurre en nuestras simulaciones. 

## 10.1. Distribución inicial de velocidades

Podemos hacer el gráfico de la distribución de Maxwell-Boltzmann para ver su aspecto en función de la temperatura:

```julia-repl
julia> mb(v,kT) = (1/(2π*kT))^(3/2)*4π*v^2*exp(-v^2/(2*kT))

julia> v = 0:0.1:5 # velocity range

julia> y = mb.(v,0.6) # the . applies to all elements of x

julia> plot(v,y,label="Maxwell-Bolztmann",xlabel="v",ylabel="frequency")
```

Un aspecto interesante de ver en las simulaciones es como evoluciona, y como 
converge, la distribución de velocidades de las partículas. En nuestros programas
tenemos tres diferentes maneras de generar la distribución *inicial* de velocidades: con normas distribuidas homogéneamente entre 0 y 1, con una distribución normal de normas, o todas nulas. Esto se puede ajustar con el parámetro `initial_velocities` de `Options`. Por ejemplo, velocidades iniciales con distribuciones distintas pueden ser generadas de la siguiente forma: 

```julia-repl
julia> sys = System(n=100_000,sides=[50,50,50])
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

Para ver la distribución de velocidades en cada caso, vamos a usar la función `velocity_distribution`, que recibe lee las velocidades grabadas en un archivo de trayectorias (la razón de esto va a quedar clara en seguida):

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
Las dos distribuciones tienen la misma energía cinética media, claro, que corresponde a la temperatura deseada. Las distribuciones, por otro lado, son muy distintas:

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

A distribución normal es la que se parece a la distribución de Maxwell-Boltzmann y es la default en las simulaciones. 

## 10.2. Las velocidades en equilibrio

Vamos a comparar las distribuciones iniciales con una distribución obtenida de 
una simulación. Para eso, vamos a repetir la simulación de Langevin, pero ahora 
en 3 dimensiones, la cual fue iniciada con velocidades nulas. Vamos llamar la atención para que la simulación salva la trayectoria de las velocidades también:
```julia-repl
julia> sys = System(n=100,sides=[50,50,50])
System{Point3D}:
 Number of particles = 100
 Box sides = [50.0, 50.0, 50.0]
 
julia> minimize!(sys)

julia> md_out = md_langevin(
           sys,
           Options(
               lambda=0.01,
               nsteps=100_000,
               velocity_file="vel.xyz",
               initial_velocities=:zero
           )
       )
```

Calculamos la distribución de velocidades de esta trayectoria con:
```julia-repl
julia> v_sim = velocity_distribution(sys,"vel.xyz");
Number of frames read = 10001
 Average velocity = [0.003014992500749613, -0.04277096828317029, -0.01629199255074569]
 Average velocity norm = 1.2234527952231138
 Average kinetic energy = 0.8819808992851449
```

Haga el gráfico de esta distribución de velocidades: 
```julia-repl
julia> plot(
           v_sim,
           label="Langevin",
           xlabel="velocity norm",
           ylabel="frequency",
           linewidth=2
       )
```

Compare, visualmente, el resultado con lo que es esperado para la distribución de Maxwell-Boltzmann.

## 10.3. Ajustando los datos 

La curva de `v_sim` puede ser ajustada para ver como bien, o mal, corresponde a la distribución de Maxwell-Boltzmann. Esta tiene la forma general

$$f(v)dv = A v^2 e^{-v^2/2kT} dv$$

Esta función tiene apenas dos parámetros ajustables: $A$ y $kT$ (o simplemente $T$, pero aquí mantenemos la constante junto con la 
temperatura). 

Podemos ajustar esta curva a los datos obtenidos usando un software de ajuste, o por ejemplo el paquete [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl):

### Ajustando el modelo

```julia-repl
julia> using LsqFit

julia> @. f(v,p) = p[1]*v^2*exp(-v^2/(2*p[2]))

julia> x = v_sim[1]; y = v_sim[2]; # to make code clearer

julia> p0 = rand(2) # punto inicial

julia> fit = curve_fit(f, x, y, p0)
```

Y vemos como bien los datos son ajustados con el modelo, calculando la suma de cuadrados de los residuos:

```julia-repl
julia> sum(fit.resid.^2)
2.016208460480086e-6
```

Y podemos visualizar el modelo ajustado gráficamente, con:
```julia-repl
julia> pars = coef(fit)

julia> yfit = [ f(xi,pars) for xi in x ];

julia> plot(
           [ (x,y), (x,yfit) ],
           label=[ "Simulation"  "Fit" ],
           xlabel="velocity norm",
           ylabel="frequency",
           linewidth=2
       )
```

```@raw html
<center>
<img src=../figures/velocities_fit.svg>
</center>
```

Note que el error es muy bajo, indicando que las velocidades de la simulación se ajustan bien a la distribución de Maxwell-Boltzmann. 
Asimismo, el segundo parámetro del ajuste tiene que ser el $kT$ que fue definido para la simulación:
```julia-repl
julia> coef(fit)[2]
 0.5852713287386224
```

## 10.4. Código completo resumido

```julia
using FundamentosDMC, Plots, LsqFit

# system setup
sys = System(n=100,sides=[50,50,50])
minimize!(sys)

# simulation
md_out = md_langevin(
    sys,
    Options(
        lambda=0.01,
        nsteps=100_000,
        velocity_file="vel.xyz",
        initial_velocities=:zero
    )
)

# plot velocity distribution
v_sim = velocity_distribution(sys,"vel.xyz");
plot(
    v_sim,
    label="Langevin",
    xlabel="velocity norm",
    ylabel="frequency",
    linewidth=2
)

# Fit model
@. f(v,p) = p[1]*v^2*exp(-v^2/(2*p[2]))
x = v_sim[1]; y = v_sim[2]; # to make code clearer
p0 = rand(2) # punto inicial
fit = curve_fit(f, x, y, p0)

# compare fit with simulation
pars = coef(fit)
yfit = [ f(xi,pars) for xi in x ]
plot(
    [ (x,y), (x,yfit) ],
    label=[ "Simulation"  "Fit" ],
    xlabel="velocity norm",
    ylabel="frequency",
    linewidth=2
)
savefig("./velocities_fit.pdf")

# Check kT
coef(fit)[2]
```










