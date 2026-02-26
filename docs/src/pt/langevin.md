# Termostato de Langevin

A função [md-langevin.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-langevin.jl) implementa o controle de temperatura de Langevin.

O termostato de Langevin consiste na suposição de que cada partícula real do sistema está imersa em um fluido de partículas muito menores, que desaceleram por fricção e, ao mesmo tempo, exercem forças aleatórias sobre as partículas reais. A desaceleração por fricção é proporcional à velocidade, e os choques aleatórios causam variações instantâneas nas velocidades das partículas reais. A trajetória de uma partícula real é, assim, propagada modificando as forças e velocidades. As forças são modificadas pela introdução da fricção,

$$\vec{F}(t) = -\nabla V[\vec{x}(t)] - \lambda \vec{v}(t)$$

E para que os choques aleatórios das partículas do banho sejam instantâneos (provoquem mudanças instantâneas nas velocidades), uma opção é introduzi-los como mudanças na velocidade:

$$v(t+\Delta t) = v(t) + a(t)\Delta t + \sqrt{2\lambda kT \Delta t}\delta(t)$$

Onde $\delta(t)$ é uma variável aleatória Gaussiana com média zero e desvio padrão 1. A relação entre o coeficiente de fricção, $\lambda$, que desacelera as partículas, e a intensidade dos choques estocásticos, $\sqrt{2\lambda kT}$ é resultado do teorema de flutuação-dissipação, que descreve o movimento Browniano.

# 6.1. Controle de parâmetros e termostatização

O coeficiente de fricção, $\lambda$, controla o comportamento de uma dinâmica de Langevin. Neste caso, o programa inicializa as velocidades em zero, para que o efeito do termostato seja destacado. Vamos imprimir a trajetória com mais frequência para notar o efeito da fricção no início:

Por exemplo,
```julia-repl
julia> out = md_langevin(
           sys,
           Options(lambda=0.1,dt=0.05,nsteps=20_000,iprintxyz=5,initial_velocities=:zero)
       );
```

Execute o programa com diversos parâmetros, em particular estes:

| $\lambda$ | $\Delta t$ |
|:---------:|:----------:|
| 0.001 |  0.05 |
| 0.01  |  0.05 |
| 10.   |  0.05 |
| 10.   |  0.001|
|       |       |

Após cada execução, observe os gráficos de energia e as trajetórias. Discuta se a temperatura chegou ao valor desejado (energia cinética igual a 60 unidades), e se a energia total é em média constante. Observe o movimento das partículas em cada trajetória.
A consistência do termostato depende do passo de integração, em particular para acoplamentos grandes.

## 6.2. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_langevin(
    sys,
    Options(lambda=0.1,dt=0.05,nsteps=20_000,iprintxyz=5,initial_velocities=:zero)
)
plot(
    out,ylim=[-100,100],
    label=["Potential" "Kinetic" "Total" "Temperature"],
    xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```
