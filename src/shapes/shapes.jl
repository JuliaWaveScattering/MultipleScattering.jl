"""
Abstract idea which defines the external boundary of object. Two objects have
the same shape if they are congruence.
"""
abstract type Shape{T<:AbstractFloat,Dim} end

"Origin of a shape, typically the center"
origin(shape::Shape) = shape.origin

# For two different shapes, the answer is false. For concrete types this
# function must be overloaded
"Returns true if two shapes are the same, ignoring their origin."
congruent(s1::Shape, s2::Shape) = false

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
"Name of a shape"
name

"The radius of a circle which contains the shape."
outer_radius

"Volume of a shape"
volume

"Returns whether an object (2nd arg) is inside a shape (1st arg)"
inside

"""
Returns Dim functions which accept a boundary coordinate (0<=t<=1)to trace outer
boundary of shape.
"""
boundary_functions
