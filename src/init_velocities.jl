"""

```
init_velocities(sys::System{T},opt::Options) where T
```

Initialize velocities with random distribution.

"""
function init_velocities(sys::System{T},opt::Options) where T
  @unpack n = sys
  @unpack kT = opt
  # type of random distribution to be created
  if opt.initial_velocities == :normal
    v = randn(T,n)
  elseif opt.initial_velocities == :flat
    v = rand(T,n)
    @. v = random_rotate(v/norm(v))
    @. v = rand()*v
  elseif opt.initial_velocities == :zero
    v = zeros(T,n)
  else
    error("Invalid option for Options().initial_velocities.")
  end
  # Remove possible remaining drift 
  remove_drift!(v)
  # Adjust average to desidred temperature  for i in 1:n
  kavg = kinetic(v)/n
  @. v = sqrt((dim(T)*kT/2)/kavg)*v
  return v
end

"""

```
function random_rotate(v)
```

Randomly rotate vector `v`, which might be a `Point2D` or a `Point3D`

"""
function random_rotate(v::Point2D)
  θ = 2π*rand()
  M = @SMatrix[ cos(θ) -sin(θ)
                sin(θ)  cos(θ) ]
  return M*v
end

function random_rotate(v::Point3D)
  θ = 2π*rand()
  ϕ = 2π*rand()
  ψ = 2π*rand()
  M = @SMatrix[ cos(ψ)*cos(θ)*cos(ϕ)-sin(ψ)*sin(ϕ)   -cos(ψ)*cos(θ)*sin(ϕ)-sin(ψ)*cos(ϕ)  cos(ψ)*sin(θ)
                sin(ψ)*cos(θ)*cos(ϕ)+cos(ψ)*sin(ϕ)   -sin(ψ)*cos(θ)*sin(ϕ)+cos(ψ)*cos(ϕ)  sin(ψ)*sin(θ)
                    -sin(θ)*cos(ϕ)                          sin(θ)*sin(ϕ)                    cos(θ)      ]
  return M*v
end
