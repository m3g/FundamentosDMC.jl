
# Fundamentals of Classical Molecular Dynamics

```@raw html
<br>
Leandro Martínez
<br>
<i>Instituto de Quimica, Universidade Estadual de Campinas (UNICAMP)</i>
<br>
<a href=http://m3g.iqm.unicamp.br>http://m3g.iqm.unicamp.br</a>
<br>
<br>
```

This tutorial contains explanations for running and analyzing
Molecular Dynamics and Monte Carlo simulations of a simple
two-dimensional system. The goal is for the student to become familiar
with various technical details involved in performing
simulations and their limitations.

## 1.1. Installing FundamentosDMC.jl

Install the [Julia](https://julialang.org) interpreter using the
[juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager) tool.
Julia version 1.9 will be installed.

!!! tip
    - By pressing `]` at the Julia prompt, you will access the package manager (a blue prompt). A sequence like `julia> ] add Package` is **not** meant to be copied and pasted as-is.

    - All other commands in this tutorial can be copied and pasted directly, including the `julia>` prompt, which is removed automatically.

Run the interpreter. We will create an environment where we will install the necessary packages:

```julia-repl
julia> ] activate @FundamentosDMC
```
(the `]` will take you to the package management prompt, "`(@v1.9) pkg>`", where `activate ...` should be run)
Press `backspace` to return to the `julia>` prompt.

Then, install the packages:

```julia-repl
julia> ] add https://github.com/m3g/FundamentosDMC.jl, Plots, LsqFit
```

This process may take a couple of minutes and will compile and install all packages in the `@FundamentosDMC` environment.

## 1.2. Installing VMD

We will use [VMD](https://www.ks.uiuc.edu/Research/vmd/) for trajectory visualization. Make sure to have it installed as well.

## 1.3. Code

`Julia` is a high-performance, dynamic programming language with a fairly simple syntax ([read here](https://www.nature.com/articles/d41586-019-02310-3)). All the code for this course is available in the repository
[https://github.com/m3g/FundamentosDMC.jl](https://github.com/m3g/FundamentosDMC.jl), and reading through it will be encouraged as we progress through the tutorial.

The code here is reasonably efficient for carrying out the tutorial, but is not optimized for maximum performance, so that it is easier to read.

## 1.4. Complete summary code

Packages can be installed from a script with:

```julia
import Pkg
Pkg.activate("FundamentosDMC", shared=true)
Pkg.add(url="https://github.com/m3g/FundamentosDMC.jl")
Pkg.add("Plots")
Pkg.add("LsqFit")
```



