"""

```
initial_point(n,opt::Options)
```

Create initial coordinates of `n` points, minimizes its energy, and returns it.

### Example
```julia-repl
julia> opt = Options(sides=Point(100,100)); # box with side 100.

julia> x0 = initial(1000,opt)

```

"""
function initial_point(n,opt::Options)
  # Simplify code by aliasing common variables
  @unpack sides, eps, sig = opt 
  u(x) = potential(x,opt)
  f!(f,x) = forces!(f,x,opt)

  # Creating random initial coordinates
  x = rand(Point,n)
  for i in 1:n
    x[i] = -sides/2 + sides .* x[i]
  end

  # Minimizing the energy of the initial point
  ubest = u(x)
  println(" Energy before minimization: ", ubest)

  dx = 1.
  fnorm = 1.
  xtrial = zeros(Point,n)
  f = zeros(Point,n)
  trial = 0
  while fnorm > 1.e-5
    trial += 1

    # Compute gradient (≡ -forces)
    f!(f,x)
    fnorm = 0.
    for v in f
      fnorm += norm2(v)
    end
    fnorm = fnorm/n

    # Compute trial point ( xtrial = x - (dU/dx)*dx = x + f*dx )
    for i in 1:n
      xtrial[i] = x[i] + f[i]*dx
    end
    utrial = u(xtrial)

    # If energy decreased, accept new point, if not, reject it and decrease dx
    if utrial < ubest
      @printf(
        "ACCEPTED U = %12.5e FNORM = %12.5e DX = %12.5e TRIAL = %10i\n", 
        utrial, fnorm, dx, trial 
      )
      for i in 1:n
        x[i] = xtrial[i]
      end
      ubest = utrial
      dx = dx * 2
    else
      dx = dx / 2
    end
  end

  println("Energy after minimization: ", u(x))
  return x
end
