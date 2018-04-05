"""
Abstract idea which defines the external boundary of object. Two objects have
the same shape if they are congruence.
"""
abstract type Shape{Dim,T<:AbstractFloat} end

volume(shape::Shape) = error("Volume function not implemented for this shape")

name(shape::Shape) = error("Volume function not implemented for this shape")


"Shape where boundary is fixed distance from a point"
type Circle{T} <: Shape{2,T}
    radius::T
end

function volume{T}(shape::Circle{T})
    return π * shape.radius^2
end

name(shape::Circle) = "Circle"


"Shape with 4 right-angles, defined by a height and width"
type Rectangle{T} <: Shape{2,T}
    width::T
    height::T
end

volume(rectangle::Rectangle) = rectangle.width*rectangle.height

name(shape::Rectangle) = "Rectangle"
