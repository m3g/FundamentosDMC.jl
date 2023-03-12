"""
    velocity_distribution(
      sys::System{T},
      velfile::String;
      vmax::Float64=2.,
      nbins::Int=100
    ) where T

Computes the distribution of velocities. If `T=Point3D` the system will be read as three-dimensional.

"""
function velocity_distribution(
    sys::System{T},
    velfile::String;
    vmax::Float64=2.0,
    nbins::Int=100
) where {T}

    # Initialize histogram
    @unpack n = sys
    Δv = vmax / nbins
    vvals = vmax * [(i - 1 / 2) * Δv for i in 1:nbins]
    vcount = zeros(nbins)
    vavg = zero(T)
    vavg_norm = 0.0
    kavg = 0.0

    # Reading velocities, computing distances, and adding to histogram
    nframes = 0
    open(velfile) do file
        while !(eof(file))
            readline(file) # number of atoms, ignore
            readline(file) # title, ignore
            nframes += 1
            # read atomic coordinates
            for i in 1:n
                atom_data = split(readline(file))
                v = T(ntuple(i -> parse(Float64, atom_data[i+1]), dim(T)))
                ibin = floor(Int, (norm(v) / vmax) / Δv) + 1
                if ibin <= nbins
                    vcount[ibin] += 1
                end
                vavg += v
                vavg_norm += norm(v)
                kavg += norm2(v) / 2
            end
        end
    end
    println("Number of frames read = ", nframes)

    # Normalizing by the number of frames and number of atoms
    for i in 1:nbins
        vcount[i] = vcount[i] / (nframes * n * Δv)
    end
    vavg = vavg / (nframes * n)
    println(" Average velocity = ", vavg)
    println(" Average velocity norm = ", vavg_norm / (nframes * n))
    println(" Average kinetic energy = ", kavg / (nframes * n))

    return vvals, vcount
end
