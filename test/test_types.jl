
abstract type pointy{T<:AbstractFloat} end

type point{T<:AbstractFloat} <: pointy{T}
    x::T
end

mutable struct foo{T<:AbstractFloat}
    points::Vector{pointy{T}}
end

# foo(ps::Vector{pointy{T}}) where T<:AbstractFloat = foo{T}(ps)
foo(ps::Vector{point{T}}) where T<:AbstractFloat = foo{T}(ps)


f(ps::Vector{pointy{T}}) where T<:AbstractFloat = [p.x for p in ps]
f(ps::Vector{point{T}}) where T<:AbstractFloat = [p.x for p in ps]

# import StaticArrays: SVector
# struct foo{Dim, FieldDim, T<:AbstractFloat}
#     x::Vector{SVector{Dim,T}}
#     y::SVector{FieldDim,Complex{T}}
# end
# x = [SVector(1.,2.)]
# y = SVector(1.+0.0im,2.+1.0im,3.-im)
