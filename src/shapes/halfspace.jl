"""
    Halfspace(normal::AbstractVector{T} [, origin::AbstractVector{T}=zeros()])

N-dimensional [`Shape`](@ref) defined by all points x = origin - λ * normal, for any λ > 0.
"""
struct Halfspace{T,Dim} <: Shape{T,Dim}
    normal::SVector{Dim,T} # outward pointing normal vector
    origin::SVector{Dim,T}
end

# Constructors
Halfspace(normal::AbstractVector{T}, origin::AbstractVector{T} = zeros(T,length(normal))) where {T} = Halfspace{T,length(normal)}(normal ./ norm(normal), origin)

Halfspace(normal::NTuple{Dim,T}, origin::AbstractVector{T} = zeros(T,Dim)) where {T,Dim} = Halfspace{T,Dim}(normal ./ norm(normal), origin)

name(shape::Halfspace) = "Halfspace"

outer_radius(hs::Halfspace{T}) where T = T(Inf)

# import Base.issubset
# function issubset(inner_rect::Rectangle{T}, outer_rect::Rectangle{T}) where T
#     all(topright(inner_rect) .<= topright(outer_rect)) &&
#     all(bottomleft(inner_rect) .>= bottomleft(outer_rect))
# end

import Base.in
function in(x::AbstractVector, hs::Halfspace)::Bool
    dot(x - hs.origin, hs.normal) < 0
end

import Base.(==)
function ==(h1::Halfspace{T}, h2::Halfspace{T}) where T
    h1.origin == h2.origin &&
    h1.normal == h2.normal
end

import Base.isequal
function isequal(h1::Halfspace{T}, h2::Halfspace{T}) where T
    isequal(h1.origin, h2.origin) &&
    isequal(h1.normal, h2.normal)
end

function iscongruent(h1::Halfspace{T}, h2::Halfspace{T}) where T
    h1.normal  == h2.normal
end

function congruent(h::Halfspace{T}, x) where T
    Rectangle{T}(h.normal, x)
end

function boundary_functions(h::Halfspace{T,2}, scale::T = 10.0) where T
    surface_vec = [-h.normal[2], h.normal[1]]

    function x(t)
        check_boundary_coord_range(t)
        h.origin[1] + t * scale * surface_vec[1]
    end

    function y(t)
        check_boundary_coord_range(t)
        h.origin[2] + t * scale * surface_vec[2]
    end

    return x, y
end
