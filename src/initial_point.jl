"""

```
initial_point(opt::Options)
```

Create initial point, minimizes its energy, and returns it.

### Example
```julia-repl
julia> opt = Options(n=1000,sides=Point(100,100)); # 1000 particles with side 100.

julia> x0 = initial(opt)

```

"""
function initial_point(opt::Options)
  @unpack n, sides, eps, sig = opt 
  u(x) = potential(x,opt)
  f!(x,f) = forces!(x,f,opt)

  # Creating random initial coordinates
  x = rand(Point,n)
  for i in 1:n
    x[i] = -sides/2 .+ sides .* x[i]
  end

  # Minimizing the energy of the initial point
  ulast = u(x)
  println(" Energy before minimization: ", ulast)

  dx = 1.
  fnorm = 1.
  xtrial = zeros(Point,n)
  f = zeros(Point,n)
  while fnorm > 1.e-5

    # Compute gradient
    f!(x,f)
    fnorm = 0.
    for i in 1:n
      fnorm = fnorm + f[i].x^2 + f[i].y^2
    end
    fnorm = fnorm / n

    # Compute trial point ( xtrial = x - (dU/dx)*dx )
    for i in 1:n
      xtrial[i] = x[i] + f[i]*dx
      xtrial[i] = image(xtrial[i],sides)
    end
    ustep = u(xtrial)

    # If energy decreased, accept, if not, reject and decrease dx
    if ustep < ulast
      for i in 1:n
        x[i] = xtrial[i]
      end
      ulast = ustep
      dx = dx * 2
      println("U = $ustep (accepted, fnorm = $fnorm)")
    else
      dx = dx / 2
      println("U = $ustep (rejected, fnorm = $fnorm)")
    end

  end

  print("Energy after minimization: ", u(x))
end
