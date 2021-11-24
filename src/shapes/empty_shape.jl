"""
    EmptyShape{Dim}

An empty shape with no points.
"""
struct EmptyShape{Dim} <: Shape{Dim} end

EmptyShape(s::Shape{Dim}) where {Dim} = EmptyShape{Dim}()

name(shape::EmptyShape) = "EmptyShape"

outer_radius(c::EmptyShape{Dim}) where {Dim} = 0
volume(shape::EmptyShape{Dim}) where {Dim} = 0

import Base.issubset
issubset(c1::Shape, c2::EmptyShape) = false
issubset(c1::EmptyShape, c2::Shape) = false

import Base.in
in(x::AbstractVector, c::EmptyShape) = false
