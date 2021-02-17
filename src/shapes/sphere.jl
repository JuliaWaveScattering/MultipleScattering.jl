"""
    Sphere([origin=zeros(),] radius)

A [`Shape`](@ref) where boundary is a fixed distance from the origin. In 2D this is a circle, in 3D the usual sphere, and in higher dimensions if difficult to visualise.
"""
struct Sphere{T,Dim} <: Shape{T,Dim}
    origin::SVector{Dim,T}
    radius::T
end

# Alternate constructors, where type is inferred naturally
Sphere(origin::NTuple{Dim}, radius::T) where {T,Dim} = Sphere{T,Dim}(origin, radius)
Sphere(origin::AbstractVector, radius::T) where {T} = Sphere{T,length(origin)}(origin, radius)

Sphere(Dim, radius::T) where {T} = Sphere{T,Dim}(zeros(T,Dim), radius)

Circle(origin::Union{AbstractVector{T},NTuple{2,T}}, radius::T) where T <: AbstractFloat = Sphere{T,2}(origin, radius::T)
Circle(radius::T) where T <: AbstractFloat = Sphere(2, radius::T)
Sphere(radius::T) where {T} = Sphere{T,3}(zeros(T,3), radius)

name(shape::Sphere) = "Sphere"
name(shape::Sphere{T,2}) where T = "Circle"

outer_radius(sphere::Sphere) = sphere.radius
volume(shape::Sphere{T,3}) where T = 4//3 * π * shape.radius^3
volume(shape::Sphere{T,2}) where T = π * shape.radius^2

# bounding_box(sphere::Sphere{T,3}; kws...) where T = bounding_box(Circle(sphere; kws...))
function bounding_box(sphere::Sphere{T,Dim}) where {T,Dim}
    return Box(origin(sphere), T(2)*sphere.radius .* ones(T,Dim))
end

import Base.in
function in(x::AbstractVector, sphere::Sphere)::Bool
    norm(origin(sphere) .- x) <= sphere.radius
end

import Base.issubset
function issubset(inner_sphere::Sphere{T,Dim}, outer_sphere::Sphere{T,Dim}) where {T,Dim}
    norm(origin(outer_sphere) - origin(inner_sphere)) <= outer_sphere.radius - inner_sphere.radius
end

function issubset(sphere::Sphere, box::Box)
    sphere_box = bounding_box(sphere)
    return issubset(sphere_box,box)
end

function issubset(box::Box, sphere::Sphere)
    issubset(Sphere(origin(box),outer_radius(box)), sphere)
end

import Base.(==)
function ==(c1::Sphere, c2::Sphere)
    c1.origin == c2.origin &&
    c1.radius == c2.radius
end

import Base.isequal
function isequal(c1::Sphere, c2::Sphere)
    isequal(c1.origin, c2.origin) &&
    isequal(c1.radius, c2.radius)
end

function iscongruent(s1::Sphere, s2::Sphere)
    s1.radius == s2.radius
end

function congruent(s::Sphere, x)
    Sphere(x, s.radius)
end

function Circle(sphere::Sphere{T}; y = sphere.origin[2]) where T
    if abs(y - sphere.origin[2]) > sphere.radius
        return EmptyShape(sphere)
    else
        r = sqrt(sphere.radius^2 - (sphere.origin[2] - y)^2)
        return Sphere{T,2}(sphere.origin[[1,3]], r)
    end
end

# boundary_functions(sphere::Sphere; kws...) = boundary_functions(Circle(sphere; kws...))



function boundary_functions(sphere::Sphere{T,3}) where T

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

function boundary_functions(circle::Sphere{T,2}) where T

    function x(t)
        check_boundary_coord_range(t)
        circle.radius * cos(2π * t) + origin(circle)[1]
    end

    function y(t)
        check_boundary_coord_range(t)
        circle.radius * sin(2π * t) + origin(circle)[2]
    end

    return x, y
end
