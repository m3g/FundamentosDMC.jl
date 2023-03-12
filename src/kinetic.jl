"""
  kinetic(v)

Function that computes the kinetic energy given a velocity vector `v`.

"""
function kinetic(v)
    kinetic = 0.0
    for p in v
        kinetic = kinetic + norm2(p) / 2
    end
    return kinetic
end

