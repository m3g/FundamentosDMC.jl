"""

```
forces!(x::Vector{Point},f::Vector{Point},pars::Parameters)
```

Subroutine that computes the force. It modifies the input `f` vector.

"""
function forces!(f,x,opt)
  @unpack eps, sig, sides = opt
  sig6 = sig^6
  sig12 = sig6^2
  eps4 = 4*eps
  n = length(x)
  for i in 1:n
    f[i] = zero(Point)
  end
  for i in 1:n-1
    for j in i+1:n
      dx = image(x[j] - x[i], sides)
      r = norm(dx)
      r3 = r^3
      r7 = r3*r3*r
      r13 = r7*r3*r3
      drdx = -dx/r
      dudx = eps4*(12*sig12/r13 - 6*sig6/r7)*drdx
      f[i] = f[i] + dudx 
      f[j] = f[j] - dudx
    end
  end
  return f
end