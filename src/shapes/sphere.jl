"""
    Sphere([origin=zeros(),] radius)

3D [`Shape`](@ref) where boundary is a fixed distance from the origin.
"""
struct Sphere{T} <: Shape{T,3}
    origin::SVector{3,T}
    radius::T
end

# Alternate constructors, where type is inferred naturally
Sphere(origin::NTuple{3}, radius::T) where {T} = Sphere{T}(origin, radius)
Sphere(origin::Vector, radius::T) where {T} = Sphere{T}(origin, radius)
# If no position is given, assume origin is at zero
Sphere(radius::T) where {T} = Sphere{T}(SVector(zero(T),zero(T),zero(T)), radius)

name(shape::Sphere) = "Sphere"

outer_radius(sphere::Sphere) = sphere.radius
volume(shape::Sphere) = 4//3 * π * shape.radius^3

import Base.issubset
function issubset(inner_sphere::Sphere{T}, outer_sphere::Sphere{T}) where T
    norm(origin(outer_sphere) - origin(inner_sphere)) <= outer_sphere.radius - inner_sphere.radius
end

import Base.in
function in(x::AbstractVector, sphere::Sphere)::Bool
    norm(origin(sphere) .- x) <= sphere.radius
end

function iscongruent(s1::Sphere{T}, s2::Sphere{T}) where T
    s1.radius == s2.radius
end

function congruent(s::Sphere{T}, x) where T
    Sphere{T}(x, s.radius)
end

function Circle(sphere::Sphere; y = sphere.origin[2])
    if abs(y - sphere.origin[2]) > sphere.radius
        return EmptyShape(sphere)
    else
        r = sqrt(sphere.radius^2 - (sphere.origin[2] - y)^2)
        return Circle(sphere.origin[[1,3]], r)
    end
end

bounding_rectangle(sphere::Sphere; kws...) = bounding_rectangle(Circle(sphere; kws...))
# boundary_functions(sphere::Sphere; kws...) = boundary_functions(Circle(sphere; kws...))

function boundary_functions(sphere::Sphere{T}) where T

    function x(t,s=T(0))
        check_boundary_coord_range(t)
        check_boundary_coord_range(s)

        sphere.radius * sin(t * π) * cos(s * 2π) + origin(sphere)[1]
    end

    function y(t,s=T(0))
        check_boundary_coord_range(t)
        check_boundary_coord_range(s)

        sphere.radius * sin(t * π) * sin(s * 2π) + origin(sphere)[2]
    end

    function z(t,s=T(0))
        check_boundary_coord_range(t)
        check_boundary_coord_range(s)

        sphere.radius * cos(t * π) + origin(sphere)[3]
    end

    return x, y, z
end
