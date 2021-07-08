# Initialize velocities
function velocities(n,opt::Options{T}) where T
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
    v[i] = sqrt(opt.kavg_target/kavg)*v[i]
  end
  return v
end