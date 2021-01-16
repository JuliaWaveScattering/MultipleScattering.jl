"""
    Circle([origin=zeros(),] radius)

2D [`Shape`](@ref) where boundary is a fixed distance from the origin.
"""
struct Circle{T} <: Shape{T,2}
    origin::SVector{2,T}
    radius::T
end

# Alternate constructors, where type is inferred naturally
Circle(origin::NTuple{2}, radius::T) where {T} = Circle{T}(origin, radius)
Circle(origin::Vector, radius::T) where {T} = Circle{T}(origin, radius)
# If no position is given, assume origin is at zero
Circle(radius::T) where {T} = Circle{T}(SVector(zero(T),zero(T)), radius)

import Base.show
function show(io::IO, c::Circle)
    # Circle template parameters can be determined entirely from the origin and radius so we do not need to print them
    # Print is the style of the first constructor
    write(io, "Circle($(c.origin.data), $(c.radius))")
    return
end

name(shape::Circle) = "Circle"

outer_radius(c::Circle) = c.radius
volume(shape::Circle) = π * shape.radius^2

import Base.issubset
function issubset(inner_circle::Circle{T}, outer_circle::Circle{T}) where T
    norm(origin(outer_circle) - origin(inner_circle)) <= outer_circle.radius - inner_circle.radius
end

function issubset(circle::Circle{T}, rect::Box{T,2}) where {T}
    
    all((origin(circle) .- circle.radius) .>= bottomleft(rect)) &&
    all((origin(circle) .+ circle.radius) .<= topright(rect))
end

function issubset(rect::Rectangle{T}, circle::Circle{T}) where {T}
    issubset(Circle(origin(rect),outer_radius(rect)), circle)
end

import Base.in
function in(x::AbstractVector, circle::Circle)::Bool
    norm(origin(circle) .- x) <= circle.radius
end

import Base.(==)
function ==(c1::Circle{T}, c2::Circle{T}) where T
    c1.origin == c2.origin &&
    c1.radius == c2.radius
end

import Base.isequal
function isequal(c1::Circle{T}, c2::Circle{T}) where T
    isequal(c1.origin, c2.origin) &&
    isequal(c1.radius, c2.radius)
end

function iscongruent(c1::Circle{T}, c2::Circle{T}) where T
    c1.radius == c2.radius
end

function congruent(c::Circle{T}, x) where T
    Circle{T}(x, c.radius)
end

function bounding_rectangle(circle::Circle)
    return Rectangle(origin(circle), 2*circle.radius, 2*circle.radius)
end

function boundary_functions(circle::Circle{T}) where T

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
