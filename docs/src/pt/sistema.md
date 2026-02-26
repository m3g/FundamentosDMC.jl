# Sistema simulado

A simulação é de um fluido de 100 partículas (monoatômicas) que interagem por um potencial de Lennard-Jones, em um sistema bidimensional, periódico.

$V = 4\varepsilon \left( \frac{\sigma^{12}}{r^{12}} - \frac{\sigma^6}{r^{6}} \right)$

Abra o arquivo [potential.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/potential.jl) e entenda a implementação do cálculo da energia potencial. Note que o cálculo depende de 3 parâmetros: $\varepsilon$, $\sigma$ e o tamanho do sistema periódico. Os parâmetros estão definidos na estrutura de dados `opt`, de entrada (veremos mais tarde como usá-la).

O arquivo [forces.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/forces.jl) contém o cálculo das forças (o gradiente do potencial), e o arquivo [kinetic.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/kinetic.jl) contém o cálculo da energia cinética. Como o sistema usa condições periódicas de contorno, as coordenadas têm que sempre ser calculadas em relação à imagem mínima. O cálculo da imagem mínima está implementado no arquivo [image.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/image.jl). É interessante entender a implementação de cada uma dessas funções, que são comuns a todos os métodos que vamos descrever.

## 2.1. Coordenadas iniciais

Para iniciar os trabalhos, abra uma seção do `Julia` e dê o comando:
```julia-repl
julia> using FundamentosDMC
```
O sistema inicial pode ser criado aleatoriamente, usando:
```julia-repl
julia> sys = System(n=100,sides=[100,100]) 
System{Point2D}:
 Number of particles = 100
 Box sides = [100.0, 100.0]
```
que gera `100` pontos em 2 dimensões, aleatórios, com coordenadas entre `[-50,-50]` e `[50,50]`, neste caso. `Point2D` é um tipo de variável que representa um ponto em duas dimensões. Mais adiante veremos que todo o código é genérico, e podemos fazer simulações em 3 dimensões apenas modificando o tipo de variável associado.

## 2.1. Parâmetros e opções da simulação

Os parâmetros das simulações são controlados na estrutura `Options`, por exemplo, para ajustar o passo de tempo, passamos o parâmetro `dt` para a estrutura. Isso pode ser feito na chamada das funções de simulação, como veremos.

```julia-repl
julia> Options(dt=0.1)
Options
  dt: Float64 0.1
  nsteps: Int64 2000
  eps: Float64 1.0
  sig: Float64 2.0
  initial_velocities: Symbol normal
  kT: Float64 0.6
  ibath: Int64 10
  iequil: Int64 10
  tau: Int64 10
  lambda: Float64 0.1
  alpha: Float64 0.1
  printxyz: Bool true
  iprintxyz: Int64 10
  iprint: Int64 20
  minprint: Bool false
  trajectory_file: String "traj.xyz"
  velocity_file: String "vel.xyz"
```

Neste caso, ajustamos o passo de tempo manualmente, e mantivemos todas as outras opções com valores padrão. Cada um desses parâmetros será discutido oportunamente. Note que definem o tamanho, campo de força ($\varepsilon$ e $\sigma$), $kT$ (temperatura) e os nomes dos arquivos de saída.

## 2.3. Minimização da energia

Em seguida, o ponto inicial será modificado usando o Método do Gradiente para minimizar a energia. O método consiste em mover as partículas segundo a aproximação de Taylor de ordem um, na direção de descida de energia:

$\vec{x}_{i+1} = \vec{x}_i - \nabla V(\vec{x}_i) \Delta x$

Se a energia no ponto $\vec{x}_{i+1}$ é menor que a energia no ponto $\vec{x}_i$, aceita-se o ponto $\vec{x}_{i+1}$ e o processo é repetido. Se não, $\Delta x$ é diminuído ($\Delta x = \Delta x / 2$), e um novo ponto $\vec{x}_{i+1}$ é calculado. Como a aproximação deve ser uma boa aproximação nas proximidades do ponto corrente ($\vec{x}_i$), um gradiente negativo garante que a função diminui para $\Delta x$ suficientemente pequeno. O processo é interrompido quando a norma do gradiente é pequena, ou quando muitos pontos foram testados. Em mecânica, $-\nabla V = \vec{F}$, então a função que calcula o gradiente é a mesma que calcula as forças na simulação. Abra o arquivo [minimize.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/minimize.jl) para discutir como se cria o ponto inicial.

Antes de executar a minimização de energia, vamos copiar o ponto inicial, para comparação:
```julia-repl
julia> x0 = copy(sys.x0)
```

Em seguida, minimizamos a energia com a função `minimize!`:
```julia-repl
julia> minimize!(sys)
Energy before minimization: 38322.72337856496
Energy after minimization: -74.15646912098042
```

!!! note
    Em `Julia` é uma convenção que funções que modificam seus argumentos terminam em `!`. Neste caso, a função vai modificar as posições, `sys.x0`, das partículas.

Se quiser ver a progressão das energias, ative a opção `minprint`, com `minimize!(sys,Options(minprint=true))`.

Podemos ver rapidamente o que aconteceu com as partículas, colocando-as em um gráfico. Primeiro, geramos um gráfico dos pontos antes da minimização:
```julia-repl
julia> using Plots

julia> scatter(Tuple.(x0))
```
Os pontos devem estar aleatoriamente distribuídos, e em particular alguns pontos devem estar muito próximos dos outros, o que gera potenciais muito repulsivos.

Em seguida, fazemos o gráfico do ponto com energia mínima obtido:
```julia-repl
julia> scatter(Tuple.(sys.x0))
```
e notará que os pontos agora têm uma nova disposição: há pontos formando aglomerados, porque o potencial de Lennard-Jones é atrativo em distâncias longas. Mas não há mais pontos muito próximos gerando repulsões muito grandes.

Este ponto inicial de energia mínima será usado em nossas simulações.

!!! note
    Os gráficos podem ser salvos como figuras com `savefig("plot.pdf")`, por exemplo.

## 2.4. Temperatura

A temperatura do sistema é um parâmetro também definido internamente no programa (pode ser modificado à vontade, mas não o faremos). A temperatura se define a partir da energia cinética média associada a cada grau de liberdade de movimento do sistema. No caso em que todos os movimentos podem ser escritos como translações, a definição é

$\frac{1}{2}kT = \left< \frac{1}{2} m v_x^2\right>$

onde a média, feita sobre $v_x$ aqui, é equivalente se feita sobre qualquer outro grau de liberdade de translação. Em um sistema tridimensional, portanto,

$\left<\frac{1}{2}m |\vec{v}|^2 \right> = 
\left<\frac{1}{2}m \left(v_x^2 + v_y^2 + v_z^2\right) \right> = 
3\left< \frac{1}{2} m v_x^2 \right> = \frac{3}{2}kT$

que é o resultado usual.

Nossas simulações são de um sistema bidimensional. Neste caso,

$\left< \frac{1}{2}m |\vec{v}|^2 \right> = 
\left< \frac{1}{2}m \left(v_x^2 + v_y^2\right)\right> =
2\left< \frac{1}{2}m v_x^2 \right> = kT$

Em todos os códigos foi escolhido que se objetiva simular o sistema à temperatura que corresponde a $kT = 0.6$ unidades. Os sistemas simulados têm 100 partículas e são bidimensionais, portanto a energia cinética média é de $100kT=60$ unidades. As velocidades iniciais serão geradas aleatoriamente no início da simulação.

## 2.5. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
x0 = copy(sys.x0)
minimize!(sys)
scatter(Tuple.(x0))
scatter(Tuple.(sys.x0))
```

Observação: o comando `Tuple.(x)` converte o vetor de vetores em um vetor de pares (tuplas), que é corretamente interpretado pelo `Plots` como uma única série de pontos.
