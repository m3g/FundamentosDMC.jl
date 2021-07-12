"""

```
init_velocities(sys::System{T},opt::Options) where T
```

Initialize velocities with random distribution.

"""
function init_velocities(sys::System{T},opt::Options;dist=:normal) where T
  @unpack n = sys
  @unpack kavg_target = opt
  # type of random distribution to be created
  if dist == :normal
    v = randn(T,n)
  elseif dist == :flat
    v = [ -0.5 .+ rand(T) for i in 1:n ]
  else
    error("`dist` must be either `:flat` nor `:normal`")
  end
  # Remove possible remaining drift 
  vdrift = mean(v)
  # Compute kinetic energy
  kavg = kinetic(v)/n
  # Adjust average to desidred temperature (Ref() is used to avoid
  # broadcasting over the *elements* of vdrift. 
  for i in 1:n
    v[i] = sqrt(kavg_target/kavg)*(v[i] - vdrift)
  end
  return v
end
