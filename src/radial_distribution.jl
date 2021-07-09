"""

```
radial_distribution(trajfile::String,dmax::Float64=10.,nbins::Int=10,T::DataType=Point2D)
```

Computes the g(r). If `ptype=Point3D` the system will be read as three-dimensional.

"""
function radial_distribution(
  sys::System{T},
  trajfile::String;
  dmax::Float64=10.,
  nbins::Int=100
) where T

  function shell_volume(r,Δr,T::DataType)
    if T == Point2D
      return π*((r+Δr)^2 - r^2)
    elseif T == Point3D
      return (4π/3)*((r+Δr)^3 - r^3)
    else
      error("T must be Point2D or Point3D")
    end
  end

  # Initialize histogram
  @unpack n, sides = sys
  x = zeros(T,n)
  density = n / prod(sides) 
  
  Δr = dmax/nbins
  r = [ (i-1/2)*Δr for i in 1:nbins ]
  g = zeros(nbins)

  # Reading coordinates, computing distances, and adding to histogram
  file = open(trajfile,"r")
  nframes = 0
  while true
    length(readline(file,keep=true)) == 0 && break # number of atoms or end-of-file, ignore or exit
    readline(file) # title, ignore

    nframes += 1 
    # read atomic coordinates
    for i in 1:n
      atom_data = split(readline(file))
      x[i] = T(ntuple(i -> parse(Float64,atom_data[i+1]), length(x[i])))
    end

    # Compute distances and add to histogram
    for i in 1:n-1
      for j in i+1:n
        dx = image(x[j] - x[i],sides)
        idist = ceil(Int,nbins*(norm(dx)/dmax))
        if idist <= nbins 
          g[idist] += 1.
        end
      end
    end

  end
  close(file)
  println("Number of frames read = ", nframes)

  # Normalizing by the number of frames and number of atoms
  # Normalizing by the number of atoms expected in an area or volume 
  for i in 1:nbins
    g[i] = g[i] / (nframes*((n-1)/2))
    g[i] = g[i] / (density*shell_volume(r[i],Δr,T))
  end

  return r, g
end