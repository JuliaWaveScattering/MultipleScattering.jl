"""
    EmptyShape{T,Dim}

An empty shape with no points.
"""
struct EmptyShape{T,Dim} <: Shape{T,Dim} end

EmptyShape(s::Shape{T,Dim}) where {T,Dim} = EmptyShape{T,Dim}()

name(shape::EmptyShape) = "EmptyShape"

outer_radius(c::EmptyShape{T,Dim}) where {T,Dim} = zero(T)
volume(shape::EmptyShape{T,Dim}) where {T,Dim} = zero(T)

import Base.issubset
issubset(c1::Shape, c2::EmptyShape) = false
issubset(c1::EmptyShape, c2::Shape) = false

import Base.in
in(x::AbstractVector, c::EmptyShape) = false
