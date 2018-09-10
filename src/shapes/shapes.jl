"""
Abstract idea which defines the external boundary of object.
"""
abstract type Shape{T<:AbstractFloat,Dim} end

"""
    origin(shape::Shape)::SVector

Origin of shape, typically the center
"""
origin(shape::Shape) = shape.origin

"""
    iscongruent(p1::Shape, p2::Shape)::Bool

True if shapes are the same but in different positions (origins), standard mathematical definition.
"""
iscongruent(s1::Shape, s2::Shape) = false # false by default, overload in specific examples

≅(s1::Shape, s2::Shape) = iscongruent(s1, s2)

"""
    congruent(s::Shape, x)::Shape

Create shape congruent to `s` but with origin at `x`
"""
function congruent end

"Generic helper function which tests if boundary coordinate is between 0 and 1"
function check_boundary_coord_range(t)
    if t < 0 || t > 1
        throw(DomainError("Boundary coordinate must be between 0 and 1"))
    end
end

# Concrete shapes
include("rectangle.jl")
include("circle.jl")
include("time_of_flight.jl")
include("time_of_flight_from_point.jl")
include("sphere.jl")

"points on the boundary of a shape"
function boundary_points(shape::Shape{T,Dim}, num_points::Int = 4; dr = zero(T)) where {Dim,T}
    x, y = boundary_functions(shape)
    v(τ) = SVector(x(τ),y(τ)) + dr * (SVector(x(τ),y(τ)) - origin(shape))
    return [ v(τ) for τ in linspace(zero(T),one(T),num_points+1)[1:end-1] ]
end

"Returns rectangle which completely encloses the shapes"
bounding_rectangle(shape1::Shape, shape2::Shape) = bounding_rectangle([shape1, shape2])

# Create a box which bounds an array of shapes
function bounding_rectangle(shapes::Vector{S}) where S<:Shape
    boxes = bounding_rectangle.(shapes)

    min_bottomleft = min.(bottomleft.(boxes)...)
    max_topright   = max.(topright.(boxes)...)

    return Rectangle(min_bottomleft, max_topright)
end

# Docstrings
"""
    name(shape::Shape)::String

Name of a shape
"""
name

"""
    outer_radius(shape::Shape{T})::T

The radius of a circle which completely contains the shape
"""
outer_radius

"""
    volume(shape::Shape{T})::T

Volume of a shape
"""
volume

"""
    volume(shape::Shape)::NTuple{Function,Dim)

Returns Tuple of Dim Functions which define outer boundary of shape when given boundary coordinate t∈[0,1]
"""
boundary_functions
