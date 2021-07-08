var documenterSearchIndex = {"docs":
[{"location":"sistema/#Sistema-simulado","page":"Sistema simulado","title":"Sistema simulado","text":"","category":"section"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"La simulación es de un fluido de 100 partículas (mono-atómicas) que interactúan por un potencial de Lennard-Jones, en un sistema bi-dimensional, periódico.","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"V = 4epsilon left( fracsigma^12r^12 - fracsigma^6r^6 right)","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Abra el archivo potential.jl y entienda la implementación del cálculo de la energía potencial. Note que el cálculo depende de 3 parámetros: epsilon, sigma, y el tamaño del sistema periódico. Los parámetros están definidos en la estructura de datos opt, de entrada (veremos más tarde como usarla). ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"El archivo forces.jl contiene el cálculo de las fuerzas (el gradiente del potencial), y el archivo kinetic.jl contiene el cálculo de la energía cinética. Como el sistema usa condiciones periódicas de contorno, las coordenadas tienen que siempre ser calculadas en relación a la imagen mínima. El cálculo de la imagen mínima está implementado en el archivo image.jl. Es interesante entender la implementación de cada una de estas funciones, que son comunes a todos los métodos que vamos a describir. ","category":"page"},{"location":"sistema/#.1.-Coordenadas-iniciales","page":"Sistema simulado","title":"2.1. Coordenadas iniciales","text":"","category":"section"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Para inciar los trabajos, abra una sección de Julia, y dé el comando:","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> using CELFI","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"El sistema inicial puede ser creado aleatoriamente, usando: ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> sys = System(n=100,sides=[100,100]) \nSystem{Point2D}:\n Number of particles = 100\n Box sides = [100.0, 100.0]\n","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"que genera 100 puntos en 2 dimensiones, aleatórios, con coordenadas entre [-50,-50] y [50,50], en este caso. Point2D es un tipo de variable que representa un punto en dos dimensiones. Mas adelante vamos a ver que todo el código es genérico, y podemos hacer simulaciones en 3 dimensiones apenas modificando el tipo de variable asociado. ","category":"page"},{"location":"sistema/#.1.-Parámetros-y-opciones-de-la-simulación","page":"Sistema simulado","title":"2.1. Parámetros y opciones de la simulación","text":"","category":"section"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Los parámetros de las simulaciones son controlados en la estructura Options, por ejemplo, para ajuster el paso de tiempo, pasamos el parámetro dt a la estructura. Esto puede ser echo en la llamada de las funciones de simulacion, como veremos.","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> Options(dt=0.1)\nOptions\n  dt: Float64 0.1\n  nsteps: Int64 2000\n  eps: Float64 1.0\n  sig: Float64 2.0\n  kavg_target: Float64 0.6\n  ibath: Int64 20\n  iequil: Int64 200\n  tau: Int64 200\n  printxyz: Bool true\n  iprintxyz: Int64 200\n  iprint: Int64 200\n  minprint: Bool false\n  trajectory_file: String \"traj.xyz\"\n","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"En este caso, ajustamos el paso de tiempo manualmente, y mantuvimos todas las otras opciones con valores default. Cada uno de estos parámetros será discutido oportunamente. Note que definen el tamaño, campo de fuerza (epsilon y sigma), energía cinética objetivo (temperatura), y los nombres de los archivos de salida. ","category":"page"},{"location":"sistema/#.3.-Minimización-de-la-energia","page":"Sistema simulado","title":"2.3. Minimización de la energia","text":"","category":"section"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"En seguida, el punto inicial va a ser modificado usando  el Método del Gradiente para minimizar la energía. El método consiste en mover las partículas según la aproximación de Taylor de orden uno, en la dirección de descenso de energía:","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"vecx_i+1 = vecx_i - nabla V(vecx_i) Delta x","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Si la energía en el punto vecx_i+1 es menor que la energía en el punto vecx_i, se acepta el punto vecx_i+1 y el proceso es repetido. Si no, Delta x es disminuido (Delta x = Delta x  2), y un nuevo punto vecx_i+1 es calculado. Como la aproximación debe ser una buena aproximación en las cercanias del punto corriente (vecx_i), un gradiente negativo garante que la función disminuye para Delta x suficientemente pequeño. El proceso es interrumpido cuando la norma del gradiente es pequeña, o cuando demasiados puntos fueron testados. En mecánica, -nabla V = vecF, entonces la función que calcula el gradiente es la misma que calcula las fuerzas en la simulación. Abra el archivo minimize.jl para discutir como se crea el punto inicial. ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Antes de ejecutar la minimización de energia, vamos a copiar el punto inicial, para comparación:","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> x0 = copy(sys.x0)","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"En seguida, minimizamos la energia con la función minimize!:","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> minimize!(sys)\nEnergy before minimization: 38322.72337856496\nEnergy after minimization: -74.15646912098042","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"note: Note\nEn Julia es una convención que las funciones que modifican sus argumentos terminan en !. En este caso, la función va a modificar las posiciones, sys.x0, de las partículas. ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Si desea ver la progression de las energias, active la opcion minprint, con minimize!(sys,Options(minprint=true)).","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Podemos ver rapidamente que ocurrió con las particulas, colocando-las en un gráfico. Primero, generamos un gráfico de los puntos antes de la minimización:","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> using Plots\n\njulia> scatter(Tuple.(x0))","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Los puntos deben estar aleatoriamente distribuídos, y en particular algunos puntos deben estar muy cercanos a los otros, lo que genera potenciales muy repulsivos.","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"En seguida, hacemos el gráfico del punto con energia mínima obtenido:","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"julia> scatter(Tuple.(sys.x0))","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"y notará que los puntos ahora tienen una nueva disposición: hay puntos formando clusteres, porque el potencial de Lennard-Jones es atractivo en distáncias largas. Pero no hay más puntos muy cercanos generando repulsiones muy grandes.","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Este punto inicial de energia mínima será usado en nuestras simulaciones. ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"note: Note\nLos gráficos puden ser salvos como figuras con savefig(\"plot.pdf\"), por ejemplo.","category":"page"},{"location":"sistema/#.4.-Temperatura","page":"Sistema simulado","title":"2.4. Temperatura","text":"","category":"section"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"La temperatura del sistema es un parámetro también definido internamente en el programa (puede ser modificado a gusto, pero no lo haremos). La temperatura se define a partir energía cinética media asociada a cada grado de libertad de movimiento del sistema. En el caso que todos los movimientos pueden ser escritos como translaciones, la definición es","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"frac12kT = left frac12 m v_x^2right","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"donde la media, hecha sobre v_x aqui, es equivalente si hecha sobre cualquier otro grado de libertad de translación. En un sistema tridimensional, por lo tanto, ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"leftfrac12m vecv^2 right = \nleftfrac12m left(v_x^2 + v_y^2 + v_z^2right) right = \n3left frac12 m v_x^2 right = frac32kT","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"que es el resultado usual.","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Nuestras simulaciones son de un sistema bi-dimensional. En este caso,","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"left frac12m vecv^2 right = \nleft frac12m left(v_x^2 + v_y^2right)right =\n2left frac12m v_x^2 right = kT","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"En los códigos de dinámica molecular, la definición de temperatura se da, así, por la definición de la energía cinética media o, en este caso, por kT. En el código de Monte-Carlo la definición de temperatura se da por la tasa de aceptación, con la misma definición. ","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"En todos los códigos fue escogido que se objetiva simular el sistema a la temperatura que corresponde a kT = 06 unidades. Los sistemas simulados tiene 100 partículas, por lo tanto la energía cinética media es 100kT=60 unidades. Las velocidades iniciales van a ser generadas aleatoriamente al princípio de la simulación. ","category":"page"},{"location":"sistema/#.5.-Código-completo-resumido","page":"Sistema simulado","title":"2.5. Código completo resumido","text":"","category":"section"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"using CELFI, Plots\nsys = System(n=100,sides=[100,100])\nx0 = copy(sys.x0)\nminimize!(sys)\nscatter(Tuple.(x0))\nscatter(Tuple.(sys.x0))","category":"page"},{"location":"sistema/","page":"Sistema simulado","title":"Sistema simulado","text":"Observación: el comando Tuple.(x) convierte el vector de vectores en un vector de pares (tuplas), que es correctamente interpretado por Plots como una única série de puntos.","category":"page"},{"location":"isokinetic/#Control-de-temperatura-isocinético","page":"Termostato isocinético","title":"Control de temperatura isocinético","text":"","category":"section"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"La función implementada en md-isokinetic.jl implementa el control de temperatura isocinético. En este método, las velocidades son escalonadas por un parámetro lambda = sqrtT_0T a intervalos regulares, para termostatizar el sistema a la temperatura T_0.  Su ejecución demanda la definición de dos parámetros adicionales: el intervalo de tiempo entre dos escalonamientos de velocidades, y el tiempo de equilibración. El tiempo de equilibración es el tiempo dentro del cual los escalonamientos son realizados. El objetivo debe ser obtener una simulación estable, con energía cinética media adecuada a la deseada (60 unidades), después de la equilibración.","category":"page"},{"location":"isokinetic/#.1.-Control-de-parámetros-y-termalización","page":"Termostato isocinético","title":"4.1. Control de parámetros y termalización","text":"","category":"section"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"El sistema es inicializado de la misma forma que antes, esto es:","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"julia> using CELFI, Plots\n\njulia> sys = System(n=100,sides=[100,100]) \n\njulia> minimize!(sys);","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"En seguida, vamos a ejecutar la simulación, ahora con termostato isocinético, por 2000 pasos, de los cuales 200 son de equilibración. El termostato es aplicado a cada ibath pasos:","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"julia> out = md_isokinetic(sys,Options(iequil=200,nsteps=2_000,ibath=1))","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"El gráfico de las energias en función del tiempo puede ser obtenido con:","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"julia> plot(out,ylim=[-100,100],\n         label=[\"Potential\" \"Kinetic\" \"Total\" \"Temperature\"],\n         xlabel=\"step\"\n       )\n","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"Debe notar-se que la energia total no es más constante dureante el periódo de equilibración. Las energia potencial y cinética deben haber convergido un poco mejor que en la simulación sin control de temperatura, aunque esta primera simulación es muy corta. ","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"La temperatura puede ser observada con:","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"julia> plot(\n         out[:,4],\n         label=[\"Temperature\"],\n         xlabel=\"step\"\n       )","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"Note que se mantiene practicamente constante e igual a la temperatura objetivo (0.60) durante la equilibración, pero después divergéncias pueden ser observadas. Si el sistema no está equilibrado, estas divergéncias pueden ser sistemáticas. ","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"Pruebe diferentes parámetros, y entienda el efecto del tiempo de equilibración y de la frequéncia de equilibración sobre la temperatura. ","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"Una buena condición para visualizar los resultados se obtiene con ibath=50  y iequil=5_000, para nsteps=20_000. ","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"julia> out = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000))","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"En estas condiciones, normalmente, no se debe observar un desvio sistemático de las energias o de la temperatura después de la equilibración. Repita los gráficos (en el prompt de Julia, la flecha para arriba accede a los comandos anteriores).","category":"page"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"La trajectoria, traj.xyz, puede ser vista con VMD, como explicado anteriormente. ","category":"page"},{"location":"isokinetic/#.2.-Código-completo-resumido","page":"Termostato isocinético","title":"4.2. Código completo resumido","text":"","category":"section"},{"location":"isokinetic/","page":"Termostato isocinético","title":"Termostato isocinético","text":"using CELFI, Plots\nsys = System(n=100,sides=[100,100])\nminimize!(sys)\nout = md_isokinetic(sys,Options(iequil=5_000,ibath=50,nsteps=20_000))\nplot(\n  out,ylim=[-100,100],\n  label=[\"Potential\" \"Kinetic\" \"Total\" \"Temperature\"],\n  xlabel=\"step\"\n)\nplot(out[:,4],label=\"Temperature\",xlabel=\"step\")","category":"page"},{"location":"simple/#Simulación-de-Dinámica-Molecular-Microcanónica","page":"Simulación microcanónica","title":"Simulación de Dinámica Molecular Microcanónica","text":"","category":"section"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Abra el archivo md-simple.jl, que contiene el algoritmo de simulación. La simulación empieza con velocidades aleatorias, ajustadas para la media termodinámica deseada de 0.6 unidades/átomo. A esta energía cinética media le llamaremos ``temperatura''. El algoritmo de integración es Velocity-Verlet, que consiste en propagar la posiciones con","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"vecx(t+Delta t) = vecx(t) + vecv(t)Delta t + frac12veca(t)Delta t^2","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"siendo veca(t)=vecF(t)m, donde vecF(t) es la fuerza en el tiempo corriente.  La fuerza en seguida es calculada en un tiempo posterior de tiempo con","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"vecF(t+Delta t) = -nabla Vleftvecx(t)right","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"y entonces las velocidades en el instante siguiente son calculadas con","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"vecv(t+Delta t) = vecv(t) +\nfrac12left\nfracvecF(t)m+fracvecF(t+Delta t)mright","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"completando el ciclo. En este ejemplo las masas son consideradas unitarias, para simplificar. La simulación es ejecutada por nsteps pasos, con paso de integración Delta t, este siendo un parámetro de entrada, dt, definido en Options. ","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"La simulación no tiene control de temperatura o de presión. Es una propagación de la trajectória según las leyes de Newton, que deberían conservar la energia. A esto se le llama una simulación \"microcanónica\", o \"NVE\" (que conserva, en princípio, el número de partículas, el volumen y la energia total).","category":"page"},{"location":"simple/#.1.-Paso-de-integración","page":"Simulación microcanónica","title":"3.1. Paso de integración","text":"","category":"section"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Para realizar una MD simple, con un paso de integración de dt=1.0, ejecute le comando:","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"julia> out = md(sys,Options(dt=0.1))\n","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"En princípio, está previso realizar 2000 pasos de integración de las equaciones  de movimimiento.  Pruebe pasos de integración entre 1.0 y 0.01.  Note que pasa con la energía. Note que pasa con la energía cinética media, la cual fue inicializada en 0.6 unidades/átomo. Discuta la elección del paso de integración, y los valores de energía cinética obtenidos. Las simulaciones que siguen van a usar un paso de integración dt = 0.05.","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Es posible controlar la frequéncia de impresión y el número de puntos salvos en el archivo de trajectoria con las opciones iprint y iprintxyz:","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"julia> out = md(sys,Options(dt=0.1,iprint=1,iprintxyz=5))","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"El número total de pasos se controla con el parámetro nsteps.","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"La variable out que sale de la simulación es una matriz con las energias y la temperatura en cada paso de la simulación. Es probable que la simulación \"explote\" con pasos de tiempo grandes. Para visualizar este processo, podemos hacer:","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"julia> using Plots\n\njulia> plot(\n         out,ylim=[-100,100],\n         label=[\"Potential\" \"Kinetic\" \"Total\" \"Temperature\"],\n         xlabel=\"step\"\n       )","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Y obtendremos un gráfico similar a:","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"<center> \n<img src=\"../figures/exploded.png\">\n</center>","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Para pasos de tiempo menores la simulación debe conseguir llegar hasta el fin. Podemos ver el resultado nuevamente, y debe ser algo similar a: ","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"<center> \n<img src=\"../figures/not_exploded.png\">\n</center>","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Observe y trate de entender las amplitudes de las oscilaciones de las energías cinética y potencial, y las amplitudes de las oscilaciones de la energía total. A que se deben cada una de las oscilaciones? Observe como estas oscilaciones dependen del paso de integración.","category":"page"},{"location":"simple/#.2.-Visualización-de-la-trajectoria","page":"Simulación microcanónica","title":"3.2. Visualización de la trajectoria","text":"","category":"section"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Por fin, abra la trayectoria usando VMD. No es necesário salir de la sección de Julia. Al apretar ; (punto y coma) aparecerá un prompt shell>, desde el cual  es posible ejecutar a VMD, si este está instalado correctamente y disponible en el path:","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"shell> vmd traj.xyz","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Dentro de VMD, elija la representación de VDW, en ","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Graphics -> Representations -> Drawing Method -> VDW","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"y dé el comando ","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"vmd> pbc set { 100. 100. 100. } -all","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"para indicar la periodicidad del sistema.  Para representar explícitamente el sistema periódico, elija +X;+Y;-X;-Y en","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Graphics -> Representations -> Periodic","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Para salir de VMD use el comando exit, y para volver al prompt de Julia desde shell>, use backspace. ","category":"page"},{"location":"simple/#.3.-Código-completo-resumido","page":"Simulación microcanónica","title":"3.3. Código completo resumido","text":"","category":"section"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"using CELFI, Plots\n\nsys = System(n=100)\nminimize!(sys)\nout = md(sys,Options(dt=0.05))\n\nplot(out,ylim=[-100,100],\n  label=[\"Potential\" \"Kinetic\" \"Total\" \"Temperature\"],\n  xlabel=\"step\"\n)","category":"page"},{"location":"simple/","page":"Simulación microcanónica","title":"Simulación microcanónica","text":"Es generado el archivo traj.xyz que puede ser visualizado en VMD. ","category":"page"},{"location":"#Fundamentos-de-Mecánica-Estadística-y-Simulaciones","page":"Início","title":"Fundamentos de Mecánica Estadística y Simulaciones","text":"","category":"section"},{"location":"","page":"Início","title":"Início","text":"<center>\nLeandro Martínez\n<br>\n<i>Instituto de Quimica, Universidad de Campinas</i>\n<br>\n<a href=http://m3g.iqm.unicamp.br>http://m3g.iqm.unicamp.br</a>\n<br>\n<br>\n</center>","category":"page"},{"location":"","page":"Início","title":"Início","text":"Este tutorial contiene las explicaciones para rodar y analizar simulaciones de Dinámica Molecular y Mote-Carlo de un sistema bi-dimensional simple. El objetivo es que el estudiante entre en contacto con diversos detalles técnicos involucrados en la realización de simulaciones y sus limitaciones. ","category":"page"},{"location":"#.1.-Instalación-de-CELFI.jl","page":"Início","title":"1.1. Instalación de CELFI.jl","text":"","category":"section"},{"location":"","page":"Início","title":"Início","text":"Instale el interpretador de Julia, de https://julialang.org/.","category":"page"},{"location":"","page":"Início","title":"Início","text":"Ejecute el interpretador, e instale el paquete deste curso, usando: ","category":"page"},{"location":"","page":"Início","title":"Início","text":"julia> ] add https://github.com/m3g/CELFI.jl","category":"page"},{"location":"","page":"Início","title":"Início","text":"(el ] te llevará al prompt the gerenciamento de paquetes, \"(@v1.6) pkg>\", donde debe ser ejecutado el add ...)","category":"page"},{"location":"","page":"Início","title":"Início","text":"Instale también el paquete Plots, que va a ser usado para graficar algunos resultados:","category":"page"},{"location":"","page":"Início","title":"Início","text":"julia> ] add Plots","category":"page"},{"location":"#.2.-Instalación-de-VMD","page":"Início","title":"1.2. Instalación de VMD","text":"","category":"section"},{"location":"","page":"Início","title":"Início","text":"Usaremos VMD para visualización de las trajectories. Certifique-se de tenerlo instalado también.","category":"page"},{"location":"#.3.-Códigos","page":"Início","title":"1.3. Códigos","text":"","category":"section"},{"location":"","page":"Início","title":"Início","text":"Julia es un lenguage de programación de alta peformance, dinámico, y con sintaxis bastante simple. Todos los códigos de este curso están disponíbles en repositório  https://github.com/m3g/CELFI.jl, y su lectura será estimulada en la medida que avanzemos en la ejecución del tutorial. ","category":"page"},{"location":"","page":"Início","title":"Início","text":"Los códigos aqui son razonablemente eficientes para al realización del tutorial, pero no están optimizados para máxima peformance, para que su lectura sea más fácil. ","category":"page"},{"location":"#.4.-Código-completo-resumido","page":"Início","title":"1.4. Código completo resumido","text":"","category":"section"},{"location":"","page":"Início","title":"Início","text":"Los paquetes pueden ser instalados desde un script, con:","category":"page"},{"location":"","page":"Início","title":"Início","text":"import Pkg\nPkg.add(\"https://github.com/m3g/CELFI.jl\")\nPkg.add(\"Plots\")","category":"page"}]
}
