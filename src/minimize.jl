"""

```
minimize!(x,opt::Options)
```

Minimizes the potential energy of `x`. Modifies `x`.

### Example
```julia-repl
julia> x = [ opt.sides .* rand(Point2D) for i in 1:100 ];

julia> minimize!(x,opt)

```

"""
function minimize!(x::Vector{T},opt::Options{T}) where T
  # Simplify code by aliasing common variables
  n = length(x)
  u(x) = potential(x,opt)
  f!(f,x) = forces!(f,x,opt)

  # Minimizing the energy of the initial point
  ubest = u(x)
  println(" Energy before minimization: ", ubest)

  dx = 1.
  fnorm = 1.
  xtrial = zeros(T,n)
  f = zeros(T,n)
  trial = 0
  while fnorm > 1.e-5 && trial < 20_000
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

  for i in 1:n
    x[i] = image(x[i],opt.sides)
  end

  println("Energy after minimization: ", u(x))
  return x
end
