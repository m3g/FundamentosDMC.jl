"""

```
forces!(x::Vector{Point},f::Vector{Point},pars::Parameters)
```

Subroutine that computes the force. It modifies the input `f` vector.

"""
function forces!(x::Vector{Point},f::Vector{Point},opt)
  @unpack n, eps, sig, sides = opt
  for i in 1:n
    f[i] = Point(0,0)
  end
  for i in 1:n-1
    for j in i+1:n
      dx = x[j] - x[i]
      dx = image(dx,sides)
      r = sqrt(dx.x^2 + dx.y^2)
      drdx = -dx/r
      fdx = 4*eps*(12*sig^12/r^13 - 6*sig^6/r^7)*drdx
      f[i] = f[i] + fdx
      f[j] = f[j] - fdx
    end
  end
  return nothing
end