"""

Function that computes the potential energy

"""
function potential(x,opt)
  @unpack sides, eps, sig = opt
  n = length(x)
  sig6 = sig^6
  sig12 = sig6^2
  eps4 = 4*eps
  u = 0.
  for i in 1:n-1
    for j in i+1:n
      dx = image(x[j] - x[i],sides)
      r2 = norm2(dx)
      r6 = r2^3
      r12 = r6^2
      u += eps4*(sig12 / r12 - sig6 / r6)
    end
  end
  return u
end
