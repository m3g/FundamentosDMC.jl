"""

Function that computes the kinetic energy.

"""
function kinetic(v)
  kinetic = 0.
  for p in v
    kinetic = kinetic + 0.5*(p[1]^2 + p[2]^2)
  end
  return kinetic
end

