"""

Structure that contains a two-dimensional point, subtype of FieldVector from StaticArrays.

"""
struct Point2D <: FieldVector{2,Float64}
  x::Float64
  y::Float64
end
norm2(p::Point2D) = (p.x^2 + p.y^2)
norm(p::Point2D) = sqrt(norm2(p))

"""

Structure that contains a three-dimensional point, subtype of FieldVector from StaticArrays.

"""
struct Point3D <: FieldVector{3,Float64}
  x::Float64
  y::Float64
  z::Float64
end
norm2(p::Point3D) = (p.x^2 + p.y^2 + p.z^2)
norm(p::Point3D) = sqrt(norm2(p))
