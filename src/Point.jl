"""

Structure that contains a two-dimensional point, subtype of FieldVector from StaticArrays.

"""
struct Point <: FieldVector{2,Float64}
  x::Float64
  y::Float64
end

norm2(p::Point) = (p.x^2 + p.y^2)
norm(p::Point) = sqrt(norm2(p))
