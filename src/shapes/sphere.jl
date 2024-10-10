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
    IsotropicHelmholtz([origin=zeros(),] radius, aperture)

A [`Shape`](@ref) which represents a 2D thin-walled isotropic Helmholtz resonator.
"""

struct IsotropicHelmholtz{T,Dim} <: AbstractSphere{Dim}
    origin::SVector{Dim,T}
    radius::T
    aperture::T
end

"""
    Helmholtz([origin=zeros(),] radius, aperture)

A [`Shape`](@ref) which represents a general 2D Helmholtz resonator.
"""

struct Helmholtz{T,Dim} <: AbstractSphere{Dim}
    origin::SVector{Dim,T}
    radius::T
    inner_radius::T
    aperture::T
    orientation::T
end

# Alternate constructors, where type is inferred naturally
Sphere(origin::NTuple{Dim}, radius::T) where {T,Dim} = Sphere{T,Dim}(origin, radius)
Sphere(origin::AbstractVector, radius::T) where {T} = Sphere{T,length(origin)}(origin, radius)

IsotropicHelmholtz(origin::NTuple{2}, radius::T, aperture::T) where {T} = IsotropicHelmholtz{T,2}(origin, radius, aperture)

Helmholtz(origin::NTuple{Dim}, radius::T, inner_radius::T, aperture::T, orientation::T) where {T,Dim} = Helmholtz{T,Dim}(origin, radius, inner_radius, aperture, orientation)

Sphere(Dim, radius::T) where {T} = Sphere{T,Dim}(zeros(T,Dim), radius)

IsotropicHelmholtz(radius::T, aperture::T) where {T} = IsotropicHelmholtz{T,2}(zeros(T, 2), radius, aperture)

Helmholtz(radius::T, inner_radius::T, aperture::T) where {T} = Helmholtz{T,2}(zeros(T, 2), radius, inner_radius, aperture, zero(T))

Helmholtz(origin::NTuple{Dim}, radius::T; inner_radius::T=zero(T), aperture::T=zero(T), orientation::T=zero(T)) where {T,Dim} = Helmholtz{T,Dim}(origin, radius, inner_radius, aperture, orientation)

function IsotropicHelmholtz(Dim::Int, radius::T;
    origin::T=zeros(T, Dim),
    aperture::T=zero(T)) where {T}

    return IsotropicHelmholtz{T,Dim}(origin, radius, aperture)
end

function Helmholtz(Dim::Int, radius::T; 
        origin::AbstractVector{T} = zeros(T,Dim), 
        aperture::T = zero(T), 
        inner_radius::T = zero(T),
        orientation::T = zero(T)) where {T} 
    
    return Helmholtz{T,Dim}(origin, radius, inner_radius, aperture, orientation)
end

Circle(origin::Union{AbstractVector{T},NTuple{2,T}}, radius::T) where T <: AbstractFloat = Sphere{T,2}(origin, radius::T)
Circle(radius::T) where T <: AbstractFloat = Sphere(2, radius::T)
Sphere(radius::T) where {T} = Sphere{T,3}(zeros(T,3), radius)

name(shape::Sphere) = "Sphere"
name(shape::Sphere{T,2}) where T = "Circle"
name(shape::IsotropicHelmholtz) = "IsotropicHelmholtz"
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

function ==(c1::IsotropicHelmholtz, c2::IsotropicHelmholtz)
    c1.origin == c2.origin && c1.radius == c2.radius && c1.aperture == c2.aperture
end

function ==(c1::Helmholtz, c2::Helmholtz)
    c1.origin == c2.origin && c1.radius == c2.radius && c1.inner_radius == c2.inner_radius && c1.aperture == c2.aperture && c1.orientation == c2.orientation
end

import Base.isequal
function isequal(c1::AbstractSphere, c2::AbstractSphere)
    false
end

function isequal(c1::Sphere, c2::Sphere)
    isequal(c1.origin, c2.origin) && isequal(c1.radius, c2.radius)
end

function isequal(c1::IsotropicHelmholtz, c2::IsotropicHelmholtz)
    isequal(c1.origin, c2.origin) && isequal(c1.radius, c2.radius) && isequal(s1.aperture, s2.aperture)
end

function isequal(c1::Helmholtz, c2::Helmholtz)
    isequal(c1.origin, c2.origin) && isequal(c1.radius, c2.radius) && isequal(c1.inner_radius, c2.inner_radius) && isequal(s1.aperture, s2.aperture) && isequal(c1.orientation, c2.orientation)
end

function iscongruent(c1::AbstractSphere, c2::AbstractSphere)
    false
end

function iscongruent(s1::Sphere, s2::Sphere)
    s1.radius == s2.radius
end

function iscongruent(s1::IsotropicHelmholtz, s2::IsotropicHelmholtz)
    s1.radius == s2.radius && s1.aperture == s2.aperture
end

function iscongruent(s1::Helmholtz, s2::Helmholtz)
    s1.radius == s2.radius && s1.aperture == s2.aperture && s1.inner_radius == s2.inner_radius && s1.orientation == s2.orientation
end

function congruent(s::Sphere, x)
    Sphere(x, s.radius)
end

function congruent(s::IsotropicHelmholtz, x)
    IsotropicHelmholtz(x, s.radius, s.aperture)
end

function congruent(s::Helmholtz, x, θ)
    Helmholtz(x, s.radius, s.inner_radius, s.aperture, θ)
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



function boundary_functions(sphere::Sphere{T,3}) where T

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

function boundary_functions(helm::Helmholtz{T,2}) where T

    # width of the Helmholtz wall
    w = helm.radius - helm.inner_radius

    length_outer = 2pi*helm.radius - helm.aperture
    length_inner = 2pi*helm.inner_radius - helm.aperture 
    lengths = [0.0, w, length_outer, w, length_inner]
    
    length_total = sum(lengths)
    cum_lengths = cumsum(lengths)

    θ = helm.orientation
    a = helm.aperture
    xo = helm.origin

    # Going to calculate the 4 points defining the aperture
    vθ = [cos(θ), sin(θ)]
    va = [-sin(θ), cos(θ)]
    
    b = sqrt(helm.inner_radius^2 - (a/2)^2)
    ps = [
        xo + b * vθ + a / 2 * va, 
        xo + (b+w) * vθ + a / 2 * va, 
        xo + (b+w) * vθ - a / 2 * va, 
        xo + (b) * vθ - a / 2 * va
    ]

    function x(t)
        check_boundary_coord_range(t)

        # coord length
        s = t * length_total

        i = findlast(s .>= cum_lengths[1:4])
        s = s - cum_lengths[i]

        if i == 1
            ps[1][1] + s * cos(θ)

        elseif i == 2
            θ1 = θ + a / (2*helm.radius)  
            xo[1] + helm.radius * cos(θ1 + s / helm.radius) 

        elseif i == 3 
            ps[3][1] - s * cos(θ)    
        else
            θ1 = θ - a / (2*helm.inner_radius) 
            xo[1] + helm.inner_radius * cos(θ1 - s / helm.inner_radius)
        end    
    end

    function y(t)
        check_boundary_coord_range(t)

        # coord length
        s = t * length_total

        i = findlast(s .>= cum_lengths[1:4])
        s = s - cum_lengths[i]

        if i == 1
            ps[1][2] + s * sin(θ)

        elseif i == 2
            θ1 = θ + a / (2*helm.radius) 
            xo[2] + helm.radius * sin(θ1 + s / helm.radius) 

        elseif i == 3 
            ps[3][2] - s * sin(θ)

        else
            θ1 = θ - a / (2*helm.inner_radius)
            xo[2] + helm.inner_radius * sin(θ1 - s / helm.inner_radius)
        end    
    end

    return x, y
end
