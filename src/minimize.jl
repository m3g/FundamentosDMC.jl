"""

```
minimize!(sys::System{T},opt::Options)
```

Minimizes the potential energy of `sys.x0`. Modifies `sys.x0`.

### Example
```julia-repl
julia> sys = System(n=100,sides=[100,100])

julia> minimize!(sys,Options(tol=1e-3))
Energy before minimization: 38322.72337856496
Energy after minimization: -74.15646912098042

```

"""
function minimize!(sys::System{T},opt::Options=Options()) where T
  # Simplify code by aliasing common variables
  @unpack n, x0, sides = sys
  x = x0
  u(x) = potential(x,sys,opt)
  f!(f,x) = forces!(f,x,sys,opt)

  # Minimizing the energy of the initial point
  ubest = u(x)
  println("Energy before minimization: ", ubest)

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
      if opt.minprint
        @printf(
          "ACCEPTED U = %12.5e FNORM = %12.5e DX = %12.5e TRIAL = %10i\n", 
          utrial, fnorm, dx, trial 
        )
      end
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
    x[i] = image(x[i],sides)
  end

  println("Energy after minimization: ", u(x))
  return nothing 
end
