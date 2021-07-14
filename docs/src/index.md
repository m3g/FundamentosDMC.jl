
# Fundamentos de Dinámica Molecular Clásica

```@raw html
<br>
Leandro Martínez
<br>
<i>Instituto de Quimica, Universidad de Campinas</i>
<br>
<a href=http://m3g.iqm.unicamp.br>http://m3g.iqm.unicamp.br</a>
<br>
<br>
```

Esta tutoría contiene las explicaciones para rodar y analizar
simulaciones de Dinámica Molecular y Mote-Carlo de un sistema
bidimensional simple. El objetivo es que el estudiante entre en
contacto con diversos detalles técnicos involucrados en la realización
de simulaciones y sus limitaciones. 

## 1.1. Instalación de FundamentosDMC.jl

Instale el interpretador de `Julia`, de [https://julialang.org/](https://julialang.org).

Ejecute el interpretador, e instale el paquete de este curso, usando: 

```julia-repl
julia> import Pkg

julia> Pkg.add(url="https://github.com/m3g/FundamentosDMC.jl")
```
o, alternativamente,
```julia-repl
julia> ] add https://github.com/m3g/FundamentosDMC.jl
```
(el `]` te llevará al prompt de gerenciamiento de paquetes, "`(@v1.6) pkg>`", donde debe ser ejecutado el `add ... `)

Instale también el paquete `Plots`, que va a ser usado para graficar algunos resultados:
```julia-repl
julia> ] add Plots
```

## 1.2. Instalación de VMD

Usaremos [VMD](https://www.ks.uiuc.edu/Research/vmd/) para visualización de las trayectorias. Certifique-se de tenerlo instalado también.

## 1.3. Códigos

`Julia` es un lenguaje de programación de alto desempeño, dinámico, y con sintaxis bastante simple ([lea aquí](https://www.nature.com/articles/d41586-019-02310-3)). Todos los códigos de este curso están disponibles en repositorio 
[https://github.com/m3g/FundamentosDMC.jl](https://github.com/m3g/FundamentosDMC.jl), y su lectura será estimulada en la medida que avancemos en la ejecución de la tutoría. 

Los códigos aquí son razonablemente eficientes para al realización de la tutoría, pero no están optimizados para máximo desempeño, para que su lectura sea más fácil. 

## 1.4. Código completo resumido

Los paquetes pueden ser instalados desde un script, con:

```julia
import Pkg
Pkg.add(url="https://github.com/m3g/FundamentosDMC.jl")
Pkg.add("Plots")
```



