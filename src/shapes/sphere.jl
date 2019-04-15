"""
    Sphere([origin=zeros(),] radius)

3D [`Shape`](@ref) where boundary is a fixed distance from the origin.
"""
struct Sphere{T} <: Shape{T,3}
    origin::SVector{3,T}
    radius::T
end

# Alternate constructors, where type is inferred naturally
Sphere(origin::NTuple{3,T}, radius::T) where {T} = Sphere{T}(origin, radius)
Sphere(origin::Vector{T}, radius::T) where {T} = Sphere{T}(origin, radius)
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
