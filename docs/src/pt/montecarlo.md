# Simulações de Monte Carlo

A função [mc.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/mc.jl) implementa o método de Monte Carlo.

As simulações de Monte Carlo têm um princípio totalmente diferente das simulações de dinâmica, mas supõe-se que amostram o mesmo conjunto de configurações se as condições termodinâmicas forem as mesmas. Aqui realizaremos uma simulação de Monte Carlo e verificaremos quão semelhante ela é em relação às simulações de dinâmica molecular.

Ao contrário de MD, MC não tem tempo. Há uma geração de posições aleatórias consecutivas, que são aceitas ou não de acordo com o critério de Metrópolis,

Se $V(\vec{x}_j) \leqslant V(\vec{x}_i)$, então $P(i\to j) = 1$

Se $V(\vec{x}_j) > V(\vec{x}_i)$, então $P(i\to j) = e^{-(V_j-V_i)/kT}$

O segundo critério é, numericamente, satisfeito comparando o resultado de $e^{-(V_j-V_i)/kT}$ com o sorteio de um número aleatório entre 0 e 1. Em nossos exemplos, $kT=0.6$.

Esse procedimento gera uma sequência de configurações, que na prática tem correlação porque as novas configurações geralmente são geradas por perturbações das configurações anteriores. As perturbações devem ser escolhidas para minimizar a correlação, ao mesmo tempo em que a taxa de aceitação seja razoável. Taxas de aceitação da ordem de 20 a 30% são consideradas ideais.

## 7.1. Execução

Execute a função `mc`. Vamos variar a magnitude das perturbações. As perturbações das posições são Gaussianas, e a magnitude de entrada é o desvio padrão. O número de passos corresponde ao número de novas estruturas, não necessariamente aceitas, geradas:

```julia-repl
julia> out = mc(sys,Options(alpha=0.05,nsteps=50_000));
```

Para um número de passos de 50.000, teste diferentes perturbações, até que uma taxa de aceitação de cerca de 30% seja obtida. (Algo próximo de $0.08\textrm{\AA}$).

Uma vez escolhida a perturbação, execute o programa com número de passos de 200.000, o que implica que aproximadamente 60.000 passos serão aceitos (para uma taxa de 30%).

Observe a evolução da energia potencial, fazendo gráficos com:

```julia-repl
plot(out,ylim=[-100,100], label="Potential", xlabel="step")
```

Neste caso não faz sentido mostrar a energia cinética, que não está definida já que as partículas não têm efetivamente velocidades. Salve o gráfico em um arquivo `pdf`, para comparação posterior, usando:
```julia-repl
julia> savefig("./mc.pdf")
```

Observe a trajetória (com os mesmos comandos de antes).

## 7.2. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = mc(sys,Options(alpha=0.05,nsteps=50_000))
plot(out,ylim=[-100,100], label="Potential", xlabel="step")
```
