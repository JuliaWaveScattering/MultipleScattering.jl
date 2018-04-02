"""
Abstract idea which defines the external boundary of object. Two objects have
the same shape if they are congruence.
"""
abstract type Shape{T,Dim} end

volume(shape::Shape) = error("Volume function not implemented for shape")

name(shape::Shape) = error("Volume function not implemented for shape")


"Shape where boundary is fixed distance from a point"
type Circle{T} <: Shape{T,2}
    radius::T
end

function volume{T}(shape::Circle{T})
    return π * shape.radius^2
end

name(shape::Circle) = "Circle"


"Shape with 4 right-angles, defined by a height and width"
type Rectangle{T} <: Shape{T,2}
    width::T
    height::T
end

function volume(rectangle::Rectangle{T}) = rectangle.width*rectangle.height

name(shape::Rectangle) = "Rectangle"
