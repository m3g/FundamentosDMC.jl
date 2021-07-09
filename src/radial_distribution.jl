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
  nbins::Int=20
) where T

  function shell_volume(ibin,nbins,dmax,T::DataType)
    Δr = dmax/nbins
    r = (ibin-1)*Δr
    if T == Point2D
      return π*((r+Δr)^2 - r^2)
    elseif T == Point3D
      return (4π/3)*((r+Δr)^3 - r^3)
    end
  end

  # Initialize histogram
  @unpack n, sides = sys
  x = zeros(T,n)
  gr = zeros(nbins)
  density = n / prod(sides) 

  # Reading coordinates, computing distances, and adding to histogram
  file = open(trajfile,"r")
  nframes = 0
  while true
    length(readline(file)) == 0 && break
    readline(file)

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
          gr[idist] = gr[idist] + 1.
        end
      end
    end

  end
  close(file)
  println("Number of frames read = ", nframes)

  # Normalizing by the number of frames and number of atoms
  # Normalizing by the number of atoms expected in an area or volume 
  for i in 1:nbins
    gr[i] = gr[i] / (nframes*(n-1))
    gr[i] = gr[i] / (density*shell_volume(i,nbins,dmax,T))
  end

  return gr
end














