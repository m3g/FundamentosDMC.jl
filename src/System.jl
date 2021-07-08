"""

Structure containing the system initial coordinates and box sides.

### Example
```julia-repl
julia> sys = System(n=100,sides=[100,100])

```

"""
struct System{T}
  n::Int
  x0::Vector{T}
  sides::T
end

function System(;
  n::Int = 100,
  x0 = nothing,
  sides = nothing,
) 

  if sides === nothing
    sides = Point2D(100,100)
  else
    if length(sides) == 2
      sides = Point2D(sides...)
    elseif length(sides) == 3
      sides = Point3D(sides...)
    else
      error("Dimension of sides must be 2 or 3.")
    end
  end

  T = typeof(sides)
  if x0 === nothing
    x0 = rand(T,n)
    for i in 1:n
      x0[i] = -sides/2 + sides .* x0[i]
    end
  else
    if eltype(x0) != T
      error(" Dimension of points of in `x0` and box sides must the same. ")
    end
    n = length(x0)
  end

  return System{T}(n,x0,sides)
end

function Base.show(io::IO,sys::System)
  println(typeof(sys),":")
  println(" Number of particles = ", length(sys.x0))
  print(" Box sides = ", sys.sides)
end


