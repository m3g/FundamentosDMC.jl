"""

```
velocities(sys::System{T},opt::Options) where T
```

Initialize velocities

"""
function velocities(sys::System{T},opt::Options) where T
  @unpack n = sys
  @unpack kavg_target = opt
  # Create random velocities
  v = rand(T,n)
  # Remove possible remaining drift 
  vdrift = zero(T)
  for i in 1:n 
    vdrift += v[i]
  end
  vdrift = vdrift / n
  for i in 1:n
    v[i] -= vdrift
  end
  # Compute kinetic energy
  kavg = kinetic(v)/n
  # Adjust average to desidred temperature
  for i in 1:n
    v[i] = sqrt(kavg_target/kavg)*v[i]
  end
  return v
end