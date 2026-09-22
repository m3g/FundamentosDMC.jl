# Fundamentos de Dinâmica Molecular Clássica

```@raw html
<br>
Leandro Martínez
<br>
<i>Instituto de Química, Universidade Estadual de Campinas (UNICAMP)</i>
<br>
<a href="http://m3g.iqm.unicamp.br">http://m3g.iqm.unicamp.br</a>
<br>
<br>
```

Este tutorial apresenta explicações para executar e analisar simulações de Dinâmica Molecular e Monte Carlo de um sistema bidimensional simples. O objetivo é que o estudante tenha contato com diversos detalhes técnicos envolvidos na realização dessas simulações e compreenda suas limitações.

## 1.1. Instalação do FundamentosDMC.jl

Instale o interpretador [Julia](https://julialang.org) utilizando a ferramenta [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager). A versão 1.9 do Julia será instalada.

!!! tip
    - Ao pressionar `]` no prompt do Julia, você acessa o gerenciador de pacotes (um prompt azul). Uma sequência como `julia> ] add Package` **não** deve ser copiada e colada.

    - Todos os outros comandos deste tutorial podem ser copiados e colados diretamente, inclusive no prompt `julia>`, que é removido automaticamente.

Execute o interpretador. Vamos criar um ambiente onde instalaremos os pacotes necessários:

```julia-repl
julia> ] activate @FundamentosDMC
```
(O `]` levará você ao prompt de gerenciamento de pacotes, "`(@v1.9) pkg>`", onde deve ser executado o comando `activate ...`)
Pressione `backspace` para voltar ao prompt `julia>`.

Em seguida, instale os pacotes:

```julia-repl
julia> ] add https://github.com/m3g/FundamentosDMC.jl Plots LsqFit
```

Esse processo pode levar alguns minutos, pois irá compilar e instalar todos os pacotes no ambiente `@FundamentosDMC`.

## 1.2. Instalação do VMD

Usaremos o [VMD](https://www.ks.uiuc.edu/Research/vmd/) para visualização das trajetórias. Certifique-se de que ele está instalado em seu sistema.

## 1.3. Códigos

`Julia` é uma linguagem de programação de alto desempenho, dinâmica e com sintaxe bastante simples ([leia aqui](https://www.nature.com/articles/d41586-019-02310-3)). Todos os códigos deste curso estão disponíveis no repositório [https://github.com/m3g/FundamentosDMC.jl](https://github.com/m3g/FundamentosDMC.jl), e sua leitura será incentivada à medida que avançarmos no tutorial.

Os códigos apresentados aqui são suficientemente eficientes para os propósitos do tutorial, mas não estão otimizados para máximo desempenho, visando facilitar a compreensão.

## 1.4. Código completo resumido

Os pacotes podem ser instalados a partir de um script:

```julia
import Pkg
Pkg.activate("FundamentosDMC", shared=true)
Pkg.add(url="https://github.com/m3g/FundamentosDMC.jl")
Pkg.add("Plots")
Pkg.add("LsqFit")
```

## 1.5. Interface gráfica

Além de executar as simulações a partir do REPL do Julia, como descrito neste tutorial, o `FundamentosDMC.jl` oferece uma interface gráfica interativa, construída sobre o [GLMakie](https://github.com/MakieOrg/Makie.jl). Ela pode ser usada para acompanhar, em tempo real, a execução das simulações nos diferentes ensembles apresentados aqui — a simulação microcanônica (NVE) e as simulações NVT com os termostatos isocinético, de Berendsen e de Langevin — além de acompanhar simulações de Monte Carlo, mostrando o movimento das partículas junto com as energias potencial, cinética e total, e a temperatura. Todos os parâmetros do sistema e da simulação podem ser definidos diretamente na interface.

Para usá-la, instale e carregue o `GLMakie` além do `FundamentosDMC`:
```julia-repl
julia> import Pkg; Pkg.add("GLMakie")

julia> using FundamentosDMC, GLMakie

julia> simulate_gui()
```
