
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

Instale el interpretador de [Julia](https://julialang.org), usando la
herramienta [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager).
La version 1.9 de Julia va a ser instalada.

!!! tip
    - Al pulsear `]` en el prompt the Julia, accederás al package manager (un prompt azul). Una sequencia como `julia> ] add Package` **no** és para copiar y colar.

    - Todos los otros comandos de esta tutoria pueden ser copiados y colados directamente, incluyendo en prompt de `julia>`, que es eliminado automaticamente.

Ejecute el interpretador. Vamos a crear un ambiente donde instalaremos los paquetes necesarios:

```julia-repl
julia> ] activate @FundamentosDMC 
```
(el `]` te llevará al prompt de gerenciamiento de paquetes, "`(@v1.9) pkg>`", donde debe ser ejecutado el `activate ... `)
Pulsee `backspace` para volver al prompt `julia>`.

En seguida, instalamos los paquetes:

```julia-repl
julia> ] add https://github.com/m3g/FundamentosDMC.jl, Plots, LsqFit
```

Este proceso puede tardar un par de minutos, y va a compilar y instalar todos los paquetes en el ambiente `@FundamentosDMC`.

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
Pkg.activate("FundamentosDMC", shared=true)
Pkg.add(url="https://github.com/m3g/FundamentosDMC.jl")
Pkg.add("Plots")
Pkg.add("LsqFit")
```



