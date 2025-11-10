"""
    Polygon(origin::SVector{Dim,T}, radius::T) where {T, Dim}

A [`Shape`](@ref) defined by a set of points on the boundary `boundary_points`. The `boundary_points` should be ordered such that consecutive points are connected by edges of the polygon. The `interior_points` are a set of points that lie within the polygon and can be used to quickly determine what is inside or not. See
@doc in(x::AbstractVector, poly::Polygon).
"""
struct Polygon{T,Dim} <: Shape{Dim}
    boundary_points::Vector{SVector{Dim,T}}
    interior_points::Vector{SVector{Dim,T}}
end

function Polygon(boundary_points::Vector{V}; interior_points::Vector{V} = [mean(boundary_points)]) where {V <: AbstractVector}
    T = eltype(boundary_points[1])
    Dim = length(boundary_points[1])
    Polygon{T,Dim}(boundary_points, interior_points)
end


name(shape::Polygon) = "Polygon"

bounding_box(poly::Polygon) = Box(poly.boundary_points)

import Base.in
function in(x::AbstractVector, poly::Polygon)::Bool
    # Ray-casting algorithm for point-in-polygon test
    n = length(poly.boundary_points)
    inside = false
    j = n
    for i in 1:n
        xi, yi = poly.boundary_points[i][1], poly.boundary_points[i][2]
        xj, yj = poly.boundary_points[j][1], poly.boundary_points[j][2]
        if ((yi > x[2]) != (yj > x[2])) && (x[1] < (xj - xi) * (x[2] - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

import Base.issubset
function issubset(poly::Polygon, box::Box)
    poly_box = bounding_box(poly)
    return issubset(poly_box,box)
end

"""
    issubset(box::Box, poly::Polygon)

Returns true if the corners of the box are contained within polygon, false otherwise.
"""
function issubset(box::Box, poly::Polygon)
    return all(c ∈ poly for c in corners(box))
end