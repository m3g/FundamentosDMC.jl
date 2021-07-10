
# Simulaciones de Monte-Carlo

En la función [mc.jl](https://github.com/m3g/FundamentosDMC.jl/blob/master/src/mc.jl) está implementado el método de Monte-Carlo.

Las simulaciones de Monte-Carlo tienen un principio totalmente distinto
de las simulaciones de dinámica, pero se supone que muestrean el mismo
conjunto de configuraciones si las condiciones termodinámicas son las
mismas. Aquí realizaremos una simulación de Monte-Carlo y verificaremos
que similaridad poseen en relación a las simulaciones de dinámica
molecular. 

Al contrario de MD, MC no tiene tiempo. Hay una generación de posiciones
aleatorias consecutivas, que son aceptadas o no de acuerdo con el
criterio de Metropolis,

Si $V(\vec{x}_j) \leqslant V(\vec{x}_i)$, entonces $P(i\to j) = 1$ 

Si $V(\vec{x}_j) > V(\vec{x}_i)$, entonces $P(i\to j) = e^{-(V_j-V_i)/kT}$

El segundo criterio es, numéricamente, satisfecho comparando el
resultado de $e^{-(V_j-V_i)/kT}$ con el sorteo de un número aleatorio
entre 0 y 1. En nuestros ejemplos, $kT=0.6$.

Este procedimiento genera una secuencia de configuraciones, que en la
práctica tiene correlación porque las nuevas configuraciones son
generalmente generadas por perturbaciones de las configuraciones
anteriores. Las perturbaciones tiene que ser escogidas para minimizar la
correlación, al mismo tiempo que la tasa de aceptación sea razonable.
Tasas de aceptación del orden de 20 a 30\% son consideradas ideales. 

## 7.1. Ejecución

Ejecute la función `mc`. Vamos a variar la magnitud de las perturbaciones. Las perturbaciones de las posiciones son Gaussianas, y la magnitud de
entrada es el desvío estándar. El número de pasos corresponde al número
de nuevas estructuras, no necesariamente aceptadas, generadas:

```julia-repl
julia> out = mc(sys,Options(alpha=0.05,nsteps=50_000));
```

Para un número de pasos de 50.000, pruebe diferentes perturbaciones,
hasta que al fin una tasa de aceptación de al rededor de 30\% sea
obtenida. (Algo próximo a $0.08\textrm{\AA}$). 

Una vez elegida la perturbación, ejecute el programa con número de pasos
de 200.000, lo que implica que aproximadamente 60.000 pasos van a ser
aceptados (para una tasa de 30\%). 

Observe la evolución de la energia potencial, haciendo gráficos con: 

```julia-repl
plot(out,ylim=[-100,100], label="Potential", xlabel="step")
```

En este caso no tiene sentido mostrar la energia cinética, que no está definida ya que las partículas no tienen efectivamente velocidades. Salve el gráfico en un archivo `pdf`, para comparación posterior, usando: 
```julia-repl
julai> savefig("./mc.pdf")
```

Observe la trayectoria (con los mismos comandos de antes). 

## 7.2. Código completo resumido

```julia
using FundamentosDMC, Plots
sys = System(n=100,sides=[100,100])
minimize!(sys)
out = mc(sys,Options(alpha=0.05,nsteps=50_000))
plot(out,ylim=[-100,100], label="Potential", xlabel="step")
```
