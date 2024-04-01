"""
All abstract spheres will have an origin and a radius 
"""

abstract type AbstractSphere{Dim} <: Shape{Dim} end

"""
    Sphere([origin=zeros(),] radius)

A [`Shape`](@ref) where boundary is a fixed distance from the origin. In 2D this is a circle, in 3D the usual sphere, and in higher dimensions if difficult to visualise.
"""

struct Sphere{T,Dim} <: AbstractSphere{Dim}
    origin::SVector{Dim,T}
    radius::T
end

"""
    SphericalHelmholtz([origin=zeros(),] radius, aperture)

A [`Shape`](@ref) which represents a 2D thin-walled isotropic Helmholtz resonator.
"""

struct SphericalHelmholtz{T,Dim} <: AbstractSphere{Dim}
    origin::SVector{Dim,T}
    radius::T
    aperture::T
end

# Alternate constructors, where type is inferred naturally
Sphere(origin::NTuple{Dim}, radius::T) where {T,Dim} = Sphere{T,Dim}(origin, radius)
Sphere(origin::AbstractVector, radius::T) where {T} = Sphere{T,length(origin)}(origin, radius)
SphericalHelmholtz(origin::NTuple{2}, radius::T, aperture::T) where {T} = SphericalHelmholtz{T,2}(origin, radius, aperture)

Sphere(Dim, radius::T) where {T} = Sphere{T,Dim}(zeros(T,Dim), radius)
SphericalHelmholtz(radius::T, aperture::T) where {T} = SphericalHelmholtz{T,2}(zeros(T,2), radius, aperture)

Circle(origin::Union{AbstractVector{T},NTuple{2,T}}, radius::T) where T <: AbstractFloat = Sphere{T,2}(origin, radius::T)
Circle(radius::T) where T <: AbstractFloat = Sphere(2, radius::T)
Sphere(radius::T) where {T} = Sphere{T,3}(zeros(T,3), radius)

name(shape::Sphere) = "Sphere"
name(shape::Sphere{T,2}) where T = "Circle"
name(shape::SphericalHelmholtz) = "SphericalHelmholtz"
Symmetry(shape::AbstractSphere{Dim}) where {Dim} = RadialSymmetry{Dim}()

outer_radius(sphere::AbstractSphere) = sphere.radius
volume(shape::AbstractSphere{3}) = 4//3 * π * shape.radius^3
volume(shape::AbstractSphere{2}) = π * shape.radius^2

function Shape(sp::Sphere{T,Dim};
        addtodimensions::T = 0.0,
        vector_translation::AbstractVector{T} = zeros(T,Dim)
    ) where {T,Dim}
    Sphere(sp.origin + vector_translation, sp.radius + addtodimensions)
end

# bounding_box(sphere::Sphere{T,3}; kws...) where T = bounding_box(Circle(sphere; kws...))
function bounding_box(sphere::AbstractSphere)
    return Box(origin(sphere), fill(2*sphere.radius, dim(sphere)))
end

import Base.in
function in(x::AbstractVector, sphere::AbstractSphere)::Bool
    norm(origin(sphere) .- x) <= sphere.radius
end

import Base.issubset
function issubset(inner_sphere::AbstractSphere, outer_sphere::AbstractSphere)
    norm(origin(outer_sphere) - origin(inner_sphere)) <= outer_sphere.radius - inner_sphere.radius
end

function issubset(sphere::AbstractSphere, box::Box)
    sphere_box = bounding_box(sphere)
    return issubset(sphere_box,box)
end

function issubset(box::Box, sphere::AbstractSphere)
    issubset(Sphere(origin(box),outer_radius(box)), sphere)
end

import Base.(==)
==(c1::AbstractSphere, c2::AbstractSphere) = false

function ==(c1::Sphere, c2::Sphere)
    c1.origin == c2.origin && c1.radius == c2.radius
end

function ==(c1::SphericalHelmholtz, c2::SphericalHelmholtz)
    c1.origin == c2.origin && c1.radius == c2.radius && c1.aperture == c2.aperture
end

import Base.isequal
function isequal(c1::AbstractSphere, c2::AbstractSphere)
    false
end

function isequal(c1::Sphere, c2::Sphere)
    isequal(c1.origin, c2.origin) && isequal(c1.radius, c2.radius)
end

function isequal(c1::SphericalHelmholtz, c2::SphericalHelmholtz)
    isequal(c1.origin, c2.origin) && isequal(c1.radius, c2.radius) && s1.aperture == s2.aperture
end

function iscongruent(c1::AbstractSphere, c2::AbstractSphere)
    false
end

function iscongruent(s1::Sphere, s2::Sphere)
    s1.radius == s2.radius
end

function iscongruent(s1::SphericalHelmholtz, s2::SphericalHelmholtz)
    s1.radius == s2.radius && s1.aperture == s2.aperture
end

function congruent(s::Sphere, x)
    Sphere(x, s.radius)
end

function congruent(s::SphericalHelmholtz, x)
    SphericalHelmholtz(x, s.radius, s.aperture)
end

function Circle(sphere::AbstractSphere; y = sphere.origin[2])
    if abs(y - sphere.origin[2]) > sphere.radius
        return EmptyShape(sphere)
    else
        r = sqrt(sphere.radius^2 - (sphere.origin[2] - y)^2)
        return Sphere{number_type(sphere),2}(sphere.origin[[1,3]], r)
    end
end

# boundary_functions(sphere::AbstractSphere; kws...) = boundary_functions(Circle(sphere; kws...))



function boundary_functions(sphere::AbstractSphere{3})

    function x(t,s=0)
        check_boundary_coord_range(t)
        check_boundary_coord_range(s)

        sphere.radius * sin(t * π) * cos(s * 2π) + origin(sphere)[1]
    end

    function y(t,s=0)
        check_boundary_coord_range(t)
        check_boundary_coord_range(s)

        sphere.radius * sin(t * π) * sin(s * 2π) + origin(sphere)[2]
    end

    function z(t,s=0)
        check_boundary_coord_range(t)
        check_boundary_coord_range(s)

        sphere.radius * cos(t * π) + origin(sphere)[3]
    end

    return x, y, z
end

function boundary_functions(circle::AbstractSphere{2})

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
