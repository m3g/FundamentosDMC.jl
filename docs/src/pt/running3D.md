# Execução em 3D

Todos os códigos aqui podem gerar simulações tridimensionais. Na verdade, os únicos códigos que têm alguma dimensão definida são os que escrevem o arquivo de saída (`printxyz`) e o que calcula a função de distribuição radial (`radial_distribution`). Todos os outros códigos são genéricos, e o que fazem depende apenas do tipo de variável de entrada, ou seja, da dimensão da caixa e dos pontos.

## 9.1. Inicializando o sistema com 3 dimensões:

Basta criar o sistema em uma caixa tridimensional:

```julia-repl
julia> sys = System(n=100,sides=[50,50,50])
System{Point3D}:
 Number of particles = 100
 Box sides = [50.0, 50.0, 50.0]
```
Note que agora os pontos são do tipo `Point3D`. Efetivamente, as coordenadas, `sys.x0`, são tridimensionais:

```julia-repl
julia> sys.x0
100-element Vector{Point3D}:
 [1.7209592667160578, 2.6345231665739135, 1.9123751433362468]
 ...
 [-2.229487380572972, -2.4612177799986634, 8.742176399685668]
 [3.9099597073208123, 7.797739759305074, 4.3056981321251655]
```

A energia pode então ser minimizada como antes:
```julia-repl
julia> minimize!(sys)
Energy before minimization: 1.8502746330219281e6
Energy after minimization: -451.7109718153855
```

## 9.2. Executando e visualizando a simulação

A simulação, qualquer que seja, pode ser executada normalmente. Os códigos executados são exatamente os mesmos. Por exemplo, uma simulação com dinâmica de Langevin pode ser executada com:

```julia-repl
julia> md_out = md_langevin(
          sys,
          Options(lambda=0.01,nsteps=20_000,trajectory_file="md.xyz")
       )
```

A simulação pode ser visualizada no `VMD`, com:
```
shell> vmd md.xyz
```

E podemos visualizar a função de distribuição radial, com

```julia-repl
julia> rmd, gmd = radial_distribution(sys,"md.xyz")

julia> plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
```

## 9.3. Código completo resumido

```julia
using FundamentosDMC, Plots

sys = System(n=100,sides=[50,50,50])
minimize!(sys)

md_out = md_langevin(
    sys,
    Options(lambda=0.01,nsteps=20_000,trajectory_file="md.xyz")
)
rmd, gmd = radial_distribution(sys,"md.xyz")
plot(rmd,gmd,xlabel="r",ylabel="g(r)",label="MD")
```
