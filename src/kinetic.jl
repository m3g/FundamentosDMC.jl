"""

Function that computes the kinetic energy.

"""
function kinetic(v)
  kinetic = 0.
  for p in v
    kinetic = kinetic + norm2(p)/2
  end
  return kinetic
end

