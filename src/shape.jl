
abstract type Shape{T,Dim} end

type Circle{T} <: Shape{T,2}
    radius::T
end
