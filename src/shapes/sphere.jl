"Shape where boundary is a fixed distance from the origin"
struct Sphere{T} <: Shape{3,T}
    origin::SVector{3,T}
    radius::T
end

# Alternate constructors, where type is inferred naturally
Sphere(origin::NTuple{3,T}, radius::T) where {T} = Sphere{T}(origin, radius)
Sphere(origin::Vector{T}, radius::T) where {T} = Sphere{T}(origin, radius)

name(shape::Sphere) = "Sphere"

outer_radius(c::Sphere) = c.radius
volume(shape::Sphere) = 4//3 * π * shape.radius^3

function inside{T}(outer_sphere::Sphere{T}, inner_sphere::Sphere{T})
    norm(origin(outer_sphere) - origin(inner_sphere)) <= outer_sphere.radius - inner_sphere.radius
end

function inside(sphere::Sphere, x::AbstractVector)
    norm(origin(sphere) .- x) <= sphere.radius
end

function congruent(c1::Sphere{T}, c2::Sphere{T}) where T
    c1.radius == c2.radius
end
