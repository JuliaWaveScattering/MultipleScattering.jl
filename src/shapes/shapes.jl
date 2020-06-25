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
    ≅(p1::Shape, p2::Shape)::Bool

True if shapes are the same but in different positions (origins), standard mathematical definition.
"""
iscongruent(s1::Shape, s2::Shape) = false # false by default, overload in specific examples

# Define synonym for iscongruent ≅, and add documentation
≅(s1::Shape, s2::Shape) = iscongruent(s1, s2)
@doc (@doc iscongruent(::Shape, ::Shape)) (≅(::Shape, ::Shape))

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
include("halfspace.jl")
include("time_of_flight.jl")
include("time_of_flight_from_point.jl")
include("sphere.jl")
include("empty_shape.jl")
"""
    points_in_shape(Shape; res = 20, xres = res, yres = res,
             exclude_region = EmptyShape(region), kws...)

returns `(x_vec, region_inds)` where `x_vec` is a vector of two-dimensional points that cover a rectangle which bounds `Shape`, and `region_inds` is an array of linear indices such that `x_vec[region_inds]` are points contained `Shape`.

"""
function points_in_shape(region::Shape{T,2};
        res::Number = 20, xres::Number = res, yres::Number = res,
        exclude_region::Shape = EmptyShape(region),
        kws...) where T

    rect = bounding_rectangle(region; kws...)

    #Size of the step in x and y direction
    x_vec_step = [rect.width / xres, rect.height / yres]
    bl = bottomleft(rect)
    x_vec = [SVector{2}(bl + x_vec_step .* [i,j]) for i=0:xres, j=0:yres][:]
    region_inds = findall(x -> !(x ∈ exclude_region) && x ∈ region, x_vec)

    return x_vec, region_inds
end

function points_in_shape(region::Shape{T,3};
        res::Number = 20, xres::Number = res, zres::Number = res,
        y::T =  region.origin[2],
        exclude_region::Shape = EmptyShape(region)) where T

    rect = bounding_rectangle(region; y = y)

    #Size of the step in x and y direction
    x_vec_step = [rect.width / xres, zero(T), rect.height / zres]
    bl = bottomleft(rect)
    bl_vec = [bl[1],zero(T),bl[2]]
    x_vec = [SVector{3}(bl_vec + x_vec_step .* [i,0,j]) for i=0:xres, j=0:zres][:]
    region_inds = findall(x -> !(x ∈ exclude_region) && x ∈ region, x_vec)

    return x_vec, region_inds
end


"points on the boundary of a shape"
function boundary_points(shape::Shape{T,Dim}, num_points::Int = 4; dr = zero(T)) where {Dim,T}
    x, y = boundary_functions(shape)
    v(τ) = SVector(x(τ),y(τ)) + dr * (SVector(x(τ),y(τ)) - origin(shape))
    return [ v(τ) for τ in LinRange(zero(T),one(T),num_points+1)[1:end-1] ]
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
function name end

"""
    outer_radius(shape::Shape{T})::T

The radius of a circle which completely contains the shape
"""
function outer_radius end

"""
    volume(shape::Shape{T})::T

Volume of a shape
"""
function volume end

"""
    volume(shape::Shape)::NTuple{Function,Dim)

Returns Tuple of Dim Functions which define outer boundary of shape when given boundary coordinate t∈[0,1]
"""
function boundary_functions end

"""
    issubset(shape1, shape2)::Bool

Returns true if shape1 is entirely contained within shape2, false otherwise (also works with particles).
"""
function issubset(s1::Shape,s2::Shape) throw(MethodError(issubset, (s1,s2))) end

"""
    in(vector, shape)::Bool

Returns true if vector is in interior of shape, false otherwise.
"""
function in(::AbstractVector,::Shape) end
