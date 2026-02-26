# Controle de temperatura de Berendsen

A função [md-berendsen.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-berendsen.jl) implementa o controle de temperatura de Berendsen. Este método também é baseado no escalonamento de velocidades, mas é mais suave. As velocidades são escalonadas por

$$\lambda = \left[  
1 + \frac{\Delta t}{\tau} \left(
\frac{T_0}{T(t)} -1
\right)
\right]^{1/2}$$

onde $\Delta t$ é o passo de integração e $\tau$ é um parâmetro que define a velocidade com que o escalonamento é realizado. O escalonamento é mais suave e mais lento.

# 5.1. Controle de parâmetros e termostatização

O parâmetro $\tau$ é ajustado com a opção `tau` de `Options`. Por exemplo:

```julia-repl
julia> sys = System(n=100,sides=[100,100])

julia> minimize!(sys)

julia> out = md_berendsen(sys,Options(tau=50,iequil=500,nsteps=20_000));
```

Teste diferentes parâmetros, com `20_000` passos de simulação. Entre eles, estes:

| $\tau$ |  $i_{\mathrm{equil}}$ |
|:------:|:---------------------:|
|   50   |  500 |
|   50   |  1500|
|  300   |  1500|
|  300   |  3000|
|        |      |

Observe os gráficos de energia resultantes, usando os mesmos comandos de antes:
```julia-repl
julia> using Plots

julia> plot(
           out,ylim=[-100,100],
           label=["Potential" "Kinetic" "Total" "Temperature"],
           xlabel="step"
       )
```

Observe a suavidade, ou não, da curva de energia total. Veja se a energia cinética se aproximou da energia média desejada ($kT=60$).

## 5.2. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_berendsen(sys,Options(tau=50,iequil=500,nsteps=20_000))
plot(
    out,ylim=[-100,100],
    label=["Potential" "Kinetic" "Total" "Temperature"],
    xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```
