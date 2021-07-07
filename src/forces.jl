"""

```
forces!(x::Vector{Point},f::Vector{Point},pars::Parameters)
```

Subroutine that computes the force. It modifies the input `f` vector.

"""
function forces!(f,x,opt)
  @unpack eps, sig, sides = opt
  n = length(x)
  for i in 1:n
    f[i] = zero(Point)
  end
  for i in 1:n-1
    for j in i+1:n
      dx = image(x[j] - x[i], sides)
      r = norm(dx)
      drdx = -dx/r
      @fastpow dudx = 4*eps*(12*sig^12/r^13 - 6*sig^6/r^7)*drdx
      f[i] = f[i] + dudx 
      f[j] = f[j] - dudx
    end
  end
  return f
end