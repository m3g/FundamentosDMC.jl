
# Distribuição de velocidades

As velocidades das partículas em um sistema molecular (tridimensional) seguem a distribuição de Maxwell-Boltzmann:

$$f(v) ~dv= \left(\frac{m}{2 \pi kT}\right)^{3/2}\, 4\pi v^2 e^{ -\frac{mv^2}{2kT}} ~ dv$$

Essencialmente, essa equação diz que a probabilidade de uma partícula ter uma velocidade no intervalo $v\pm dv$ é proporcional a $v^2 e^{-v^2/kT}$. Isso se deve a dois motivos: primeiro, a probabilidade de uma partícula ter uma energia $E$ é dada pela distribuição de Boltzmann, $e^{-E/kT}$. Como a energia cinética é $mv^2/2$, o quadrado da velocidade aparece no expoente. Segundo, a *degeneração* aumenta com $v^2$ com o aumento da velocidade, porque o volume do espaço que contém os vetores de comprimento $v$ aumenta com $v^2$ em três dimensões. Isso gera o termo pré-exponencial proporcional a $v^2$.

Vamos ver se isso efetivamente ocorre em nossas simulações.

## 10.1. Distribuição inicial de velocidades

Podemos fazer o gráfico da distribuição de Maxwell-Boltzmann para ver seu aspecto em função da temperatura:

```julia-repl
julia> mb(v,kT) = (1/(2π*kT))^(3/2)*4π*v^2*exp(-v^2/(2*kT))

julia> v = 0:0.1:5 # intervalo de velocidades

julia> y = mb.(v,0.6) # o . aplica a todos os elementos de x

julia> plot(v,y,label="Maxwell-Bolztmann",xlabel="v",ylabel="frequência")
```

Um aspecto interessante de observar nas simulações é como evolui, e como
converge, a distribuição de velocidades das partículas. Em nossos programas
temos três maneiras diferentes de gerar a distribuição *inicial* de velocidades: com normas distribuídas homogeneamente entre 0 e 1, com uma distribuição normal de normas, ou todas nulas. Isso pode ser ajustado com o parâmetro `initial_velocities` de `Options`. Por exemplo, velocidades iniciais com distribuições distintas podem ser geradas da seguinte forma:

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

Para ver a distribuição de velocidades em cada caso, vamos usar a função `velocity_distribution`, que lê as velocidades gravadas em um arquivo de trajetórias (a razão disso ficará clara em breve):

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

Note que geramos `100_000` pontos, para obter uma distribuição bem amostrada.
As duas distribuições têm a mesma energia cinética média, claro, que corresponde à temperatura desejada. As distribuições, por outro lado, são bem diferentes:

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
<img src=/figures/velocities.svg>
</center>
```

A distribuição normal é a que se assemelha à distribuição de Maxwell-Boltzmann e é a padrão nas simulações.

## 10.2. As velocidades em equilíbrio

Vamos comparar as distribuições iniciais com uma distribuição obtida de
uma simulação. Para isso, vamos repetir a simulação de Langevin, agora
em 3 dimensões, iniciando com velocidades nulas. É importante garantir que a simulação salve a trajetória das velocidades também:
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

julia> @. f(v,p) = p[1]*v^2*exp(-v^2/(2*p[2]))

julia> x = v_sim[1]; y = v_sim[2]; # para clareza

julia> p0 = rand(2) # ponto inicial

julia> fit = curve_fit(f, x, y, p0)
```

E vemos o quão bem os dados são ajustados pelo modelo, calculando a soma dos quadrados dos resíduos:

```julia-repl
julia> sum(fit.resid.^2)
2.016208460480086e-6
```

E podemos visualizar o modelo ajustado graficamente, com:
```julia-repl
julia> pars = coef(fit)

julia> yfit = [ f(xi,pars) for xi in x ];

julia> plot(
           [ (x,y), (x,yfit) ],
           label=[ "Simulação"  "Ajuste" ],
           xlabel="norma da velocidade",
           ylabel="frequência",
           linewidth=2
       )
```

```@raw html
<center>
<img src=/figures/velocities_fit.svg>
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
