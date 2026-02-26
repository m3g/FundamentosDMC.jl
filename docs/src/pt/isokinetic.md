# Controle de temperatura isocinético

A função implementada em [md-isokinetic.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/md-isokinetic.jl) implementa o controle de temperatura isocinético. Neste método, as velocidades são escalonadas por um parâmetro $\lambda = \sqrt{T_0/T}$ em intervalos regulares, para termostatizar o sistema à temperatura $T_0$.

Sua execução demanda a definição de dois parâmetros adicionais: o intervalo de tempo entre dois escalonamentos de velocidades, e o tempo de equilibration. O tempo de equilibration é o tempo durante o qual os escalonamentos são realizados. O objetivo deve ser obter uma simulação estável, com energia cinética média adequada à desejada (60 unidades), após a equilibration.

# 4.1. Controle de parâmetros e termostatização

O sistema é inicializado da mesma forma que antes, ou seja:
```julia-repl
julia> using FundamentosDMC, Plots

julia> sys = System(n=100,sides=[100,100]) 

julia> minimize!(sys);
```

Em seguida, vamos executar a simulação, agora com termostato isocinético, por `2000` passos, dos quais `200` são de equilibration. O termostato é aplicado a cada `ibath` passos:

```julia-repl
julia> out = md_isokinetic(sys,Options(iequil=200,nsteps=2_000,ibath=1));
```

O gráfico das energias em função do tempo pode ser obtido com:
```julia-repl
julia> plot(out,ylim=[-100,100],
           label=["Potential" "Kinetic" "Total" "Temperature"],
           xlabel="step"
       )
```

Note que a energia total não é mais constante durante o período de equilibration. As energias potencial e cinética devem ter convergido um pouco melhor do que na simulação sem controle de temperatura, embora esta primeira simulação seja muito curta.

A temperatura pode ser observada com:
```julia-repl
julia> plot(
           out[:,4],
           label=["Temperature"],
           xlabel="step"
       )
```
Note que se mantém praticamente constante e igual à temperatura alvo (0.60) durante a equilibration, mas depois divergências podem ser observadas. Se o sistema não estiver equilibrado, essas divergências podem ser sistemáticas.

Teste diferentes parâmetros e entenda o efeito do tempo de equilibration e da frequência de equilibration sobre a temperatura.

Uma boa condição para visualizar os resultados é obtida com `ibath=50` e `iequil=5_000`, para `nsteps=20_000`.
```julia-repl
julia> out = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000));
```

Nessas condições, normalmente, não se deve observar um desvio sistemático das energias ou da temperatura após a equilibration. Repita os gráficos (no prompt do `Julia`, a seta para cima acessa os comandos anteriores).

A trajetória, `traj.xyz`, pode ser vista com `VMD`, como explicado anteriormente.

## 4.2. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000))
plot(
    out,ylim=[-100,100],
    label=["Potential" "Kinetic" "Total" "Temperature"],
    xlabel="step"
)
plot(out[:,4],label="Temperature",xlabel="step")
```
