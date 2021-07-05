"""

Structure that contains a two-dimensional point, subtype of FieldVector from StaticArrays.

"""
struct Point <: FieldVector{2,Float64}
  x::Float64
  y::Float64
end