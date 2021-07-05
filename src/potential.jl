"""

Function that computes the potential energy

"""
function potential(x,opt)
  @unpack n, sides, eps, sig = opt
  u = 0.
  for i in 1:n-1
    for j in i+1:n
      dx = image(x[j] - x[i],sides)
      r = sqrt(dx.x^2 + dx.y^2)
      u += 4*eps*(sig^12 / r^12 - sig^6 / r^6)
    end
  end
  return u
end
