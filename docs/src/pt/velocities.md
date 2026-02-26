
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
```julia-repl
```
Haga el gráfico de esta distribución de velocidades: 
julia> plot(
           label="Langevin",

temperatura). 

julia> d_flat = velocity_distribution(sys,"vflat.xyz");
Number of frames read = 1
 Average velocity = [-4.9000000006849124e-8, -1.2999999991358812e-8]
 Average velocity norm = 0.9506616041200406
 Average kinetic energy = 0.5999232099558015
```

Note que geramos `100_000` pontos, para obter uma distribuição bem amostrada.
As duas distribuições têm a mesma energia cinética média, que corresponde à temperatura desejada. As distribuições, por outro lado, são bem diferentes:

```julia-repl
julia> plot(
           [d_normal,d_flat],
           label=["Normal" "Uniforme"],
           xlabel="norma da velocidade",
           ylabel="frequência",
           linewidth=2
       )
```
que resulta em:
```@raw html
<center>
<img src=../figures/velocities.svg>
</center>
```

A distribuição normal é a que se assemelha à distribuição de Maxwell-Boltzmann e é a padrão nas simulações.

## 10.2. Velocidades em equilíbrio

Vamos comparar as distribuições iniciais com uma distribuição obtida de uma simulação. Para isso, vamos repetir a simulação de Langevin, agora em 3 dimensões, iniciando com velocidades nulas. É importante garantir que a simulação salve a trajetória das velocidades também:
```julia-repl
julia> sys = System(n=100,sides=[50,50,50])
System{Point3D}:
 Number of particles = 100
 Box sides = [50.0, 50.0, 50.0]
 
# plot velocity distribution

v_sim = velocity_distribution(sys,"vel.xyz");
           sys,
           Options(
               lambda=0.01,
               nsteps=100_000,
               velocity_file="vel.xyz",
               initial_velocities=:zero
           )
       )
```

Calculamos a distribuição de velocidades dessa trajetória com:
```julia-repl
julia> v_sim = velocity_distribution(sys,"vel.xyz");
Number of frames read = 10001
 Average velocity = [0.003014992500749613, -0.04277096828317029, -0.01629199255074569]
 Average velocity norm = 1.2234527952231138
 Average kinetic energy = 0.8819808992851449
```

Faça o gráfico dessa distribuição de velocidades:
```julia-repl
julia> plot(
           v_sim,
           label="Langevin",
           xlabel="norma da velocidade",
           ylabel="frequência",
           linewidth=2
       )
```

Compare visualmente o resultado com o que é esperado para a distribuição de Maxwell-Boltzmann.

## 10.3. Ajustando os dados

A curva de `v_sim` pode ser ajustada para verificar o quão bem ela corresponde à distribuição de Maxwell-Boltzmann. Esta tem a forma geral:

$$f(v)dv = A v^2 e^{-v^2/2kT} dv$$

Essa função tem apenas dois parâmetros ajustáveis: $A$ e $kT$ (ou simplesmente $T$, mas aqui mantemos a constante junto com a temperatura).

Podemos ajustar essa curva aos dados obtidos usando um software de ajuste, como o pacote [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl):

### Ajustando o modelo

```julia-repl
julia> using LsqFit

plot(

    v_sim,

    label="Langevin",

    xlabel="velocity norm",
```

E vemos o quão bem os dados são ajustados pelo modelo, calculando a soma dos quadrados dos resíduos:

```julia-repl
julia> sum(fit.resid.^2)
2.016208460480086e-6
```

E podemos visualizar o modelo ajustado graficamente, com:
```julia-repl
julia> pars = coef(fit)

    ylabel="frequency",

    linewidth=2
           [ (x,y), (x,yfit) ],
           label=[ "Simulação"  "Ajuste" ],
           xlabel="norma da velocidade",
           ylabel="frequência",
           linewidth=2
       )
```

```@raw html
<center>
<img src=../figures/velocities_fit.svg>
</center>
```

Note que o erro é muito baixo, indicando que as velocidades da simulação se ajustam bem à distribuição de Maxwell-Boltzmann.
Além disso, o segundo parâmetro do ajuste deve ser o $kT$ que foi definido para a simulação:
```julia-repl
julia> coef(fit)[2]
 0.5852713287386224
```

## 10.4. Código completo resumido

```julia
using FundamentosDMC, Plots, LsqFit

# configuração do sistema
sys = System(n=100,sides=[50,50,50])
minimize!(sys)

# simulação
md_out = md_langevin(
    sys,
    Options(
        lambda=0.01,
        nsteps=100_000,
        velocity_file="vel.xyz",
        initial_velocities=:zero
    )
)

# plotar distribuição de velocidades
v_sim = velocity_distribution(sys,"vel.xyz");
plot(
    v_sim,
    label="Langevin",
    xlabel="norma da velocidade",
    ylabel="frequência",
    linewidth=2
)

# Ajuste do modelo
@. f(v,p) = p[1]*v^2*exp(-v^2/(2*p[2]))
x = v_sim[1]; y = v_sim[2]; # para clareza
p0 = rand(2) # ponto inicial
fit = curve_fit(f, x, y, p0)

# comparar ajuste com simulação
pars = coef(fit)
yfit = [ f(xi,pars) for xi in x ]
plot(
    [ (x,y), (x,yfit) ],
    label=[ "Simulação"  "Ajuste" ],
    xlabel="norma da velocidade",
    ylabel="frequência",
    linewidth=2
)
savefig("./velocities_fit.pdf")

# Verificar kT
coef(fit)[2]
```
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










