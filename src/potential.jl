"""

Function that computes the potential energy

"""
function potential(x::Vector{T},sys::System{T},opt::Options) where T
  @unpack n, sides = sys
  @unpack eps, sig = opt
  u = 0.
  for i in 1:n-1
    for j in i+1:n
      dx = image(x[j] - x[i],sides)
      r = norm(dx)
      @fastpow u += 4*eps*(sig^12 / r^12 - sig^6 / r^6)
    end
  end
  return u
end
