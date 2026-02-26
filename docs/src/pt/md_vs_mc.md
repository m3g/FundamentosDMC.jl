# Comparação entre MC e MD

Vamos comparar o conjunto de estruturas geradas por uma dinâmica molecular com um conjunto de estruturas geradas por uma simulação de Monte Carlo.

## 8.1. Gerando uma boa amostragem por MD

Execute novamente o programa `md_langevin` com os parâmetros $\lambda=0.01$, por 100 mil passos, e salve a trajetória em um arquivo chamado `md.xyz`:

```julia-repl
julia> md_out = md_langevin(
           sys,
           Options(lambda=0.01,nsteps=100_000,trajectory_file="md.xyz")
       );
```

## 8.2. Gerando uma boa amostragem de MC

Execute uma longa simulação de Monte Carlo (`300_000` passos):

```julia-repl
julia> mc_out = mc(
           sys,
           Options(alpha=0.05,nsteps=300_000,trajectory_file="mc.xyz")
       );
```

## 8.3. Comparando as energias potenciais

A primeira coluna de `md_out` tem a energia potencial na simulação de dinâmica molecular. Vamos fazer o gráfico dessa energia:

```julia-repl
julia> plot(md_out[:,1],label="Potential energy- MD",xlabel="step")
```

E vamos adicionar ao mesmo gráfico a energia potencial obtida na simulação de Monte Carlo, que está em `mc_out` (note o `!` ao final de `plot!`, indicando que o gráfico anterior será modificado):

```julia-repl
julia> plot!(mc_out,label="Potential energy- MC",xlabel="step")
```

Note a similaridade, ou diferença entre os dois gráficos. Lembre-se que em nenhuma dessas simulações controlamos explicitamente a energia potencial.

## 8.4. Função de distribuição radial

Vamos comparar a estrutura média obtida usando MD com a estrutura média obtida com MC. Para isso vamos usar a função $g(r)$, chamada de "função de distribuição radial". Essa função mede a probabilidade de encontrar uma partícula a uma distância $r$ de outra partícula no sistema real, em relação a essa mesma probabilidade se não houvesse nenhuma interação entre as partículas.

No nosso caso bidimensional, o número de partículas por unidade de volume é $\rho=n/A$, onde $n=100$ é o número de partículas e $A=100^2$ é a área total do sistema simulado. O número de partículas esperado em um intervalo de distâncias entre $r$ e $r+\Delta r$ de cada partícula é, portanto, $n(r)=\rho A(r)$, onde $A(r)$ é a área de uma casca circular de raio menor $r$ e raio maior $r+\Delta r$:
```@raw html
<center>
<img src=../figures/area.svg>
</center>
```
Vemos que $A(r)=\pi (r+\Delta r)^2 - \pi r^2 \approx 2\pi r\Delta r$.
Assim, o número de partículas esperado, em média, seria $n(r)=2\pi r\Delta r\rho$, se não houvesse interações.

As interações fazem com que o número de partículas em cada distância seja diferente de uma distribuição homogênea. Se há interações favoráveis, por exemplo, a probabilidade de encontrar duas partículas próximas é maior. Essa distribuição de partículas é um dos parâmetros estruturais mais importantes.

## 8.5. Cálculo de $g(r)$

O programa [radial_distribution.jl](https://github.com/m3g/FUndamentosDMC.jl/blob/master/src/radial_distribution.jl) calcula, a partir de uma trajetória, a função $g(r)=n'(r)/n(r)$, onde $n(r)$ está definido anteriormente, e $n'(r)$ é o número médio de partículas efetivamente observado entre $r$ e $r+\Delta r$ na simulação.

Para calcular a função de distribuição radial da simulação de dinâmica molecular, execute:
```julia-repl
julia> rmd, gmd = radial_distribution(sys,"md.xyz")
```

`rmc` e `rmd` são as distâncias e as densidades relativas de partículas a cada distância, em relação à densidade do sistema.

Em seguida, obtenha o $g(r)$ da simulação de Monte Carlo:
```julia-repl
julia> rmc, gmc = radial_distribution(sys,"mc.xyz")
```

Crie um gráfico da função de distribuição da simulação de dinâmica molecular, com:
```julia-repl
julia> plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
```

Note o aumento da densidade local em distâncias curtas, que resulta das interações favoráveis entre as partículas. Note também que a densidade relativa tende a `1.0` em distâncias grandes, quando a correlação da posição entre partículas não é mais importante.

Adicione ao mesmo gráfico a função de distribuição obtida a partir da simulação de Monte Carlo:
```julia-repl
julia> plot!(rmc,gmc,xlabel="r",ylabel="g(r)",label="MC")
```
Compare. As simulações, com suas naturezas totalmente distintas, amostraram as mesmas estruturas?

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
    Options(alpha=0.05,nsteps=300_000,trajectory_file="mc.xyz")
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
