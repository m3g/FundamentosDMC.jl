"""
    forces!(x::Vector{T},f::Vector{T},opt::Options) where T

Subroutine that computes the force. It modifies the input `f` vector.

"""
function forces!(f::Vector{T}, x::Vector{T}, sys::System{T}, opt::Options) where {T}
    @unpack n, sides = sys
    @unpack eps, sig = opt
    @. f = zero(T)
    for i in 1:n-1
        for j in i+1:n
            dx = image(x[j] - x[i], sides)
            r = norm(dx)
            drdx = -dx / r
            @fastpow dudx = 4 * eps * (12 * sig^12 / r^13 - 6 * sig^6 / r^7) * drdx
            f[i] = f[i] + dudx
            f[j] = f[j] - dudx
        end
    end
    return f
end
