"""
    Plate(normal::AbstractVector{T}, width::T, [, origin::AbstractVector{T}=zeros()])

A plate defined by all the points ``\\mathbf x`` that satify ``|(\\mathbf x - \\mathbf o) \\cdot \\mathbf n| < w /2`` where ``\\mathbf n`` is the unit normal, ``\\mathbf o`` is the origin, and ``w`` is the width.
"""
struct Plate{T,Dim} <: Shape{T,Dim}
    normal::SVector{Dim,T} # a unit vector which is orthogonal to the plate
    width::T # the width
    origin::SVector{Dim,T}
end

# Constructors
Plate(normal::AbstractVector{T}, width::T, origin::AbstractVector{T} = zeros(T,length(normal))) where {T} = Plate{T,length(normal)}(normal ./ norm(normal), width, origin)

Plate(normal::NTuple{Dim,T}, width::T, origin::AbstractVector{T} = zeros(T,Dim)) where {T,Dim} = Plate{T,Dim}(normal ./ norm(normal), width, origin)

name(shape::Plate) = "Plate"
function Symmetry(shape::Plate{T,Dim}) where {T,Dim}
    return (abs(dot(shape.normal,azimuthalnormal(Dim))) == one(T)) ? PlanarAzimuthalSymmetry{Dim}() : PlanarSymmetry{Dim}()
end

volume(shape::Plate) = Inf
outer_radius(hs::Plate{T}) where T = T(Inf)

# import Base.issubset
# function issubset(inner_rect::Plate{T}, outer_rect::Plate{T}) where T
#     all(topright(inner_rect) .<= topright(outer_rect)) &&
#     all(bottomleft(inner_rect) .>= bottomleft(outer_rect))
# end

import Base.in
function in(x::AbstractVector{T}, p::Plate)::Bool where T
    abs(dot(x - p.origin, p.normal)) < p.width / T(2)
end

import Base.(==)
function ==(h1::Plate{T}, h2::Plate{T}) where T
    h1.origin == h2.origin &&
    h1.width == h2.width &&
    h1.normal == h2.normal
end

import Base.isequal
function isequal(h1::Plate{T}, h2::Plate{T}) where T
    isequal(h1.origin, h2.origin) &&
    isequal(h1.width, h2.width) &&
    isequal(h1.normal, h2.normal)
end

function iscongruent(h1::Plate{T}, h2::Plate{T}) where T
    (h1.normal  == h2.normal) && (h1.width  == h2.width)
end

function congruent(h::Plate{T}, x) where T
    Plate(h.normal, h.width, x)
end

# function boundary_functions(h::Plate{T,2}, scale::T = 10.0) where T
#     surface_vec = [-h.normal[2], h.normal[1]]
#
#     function x(t)
#         check_boundary_coord_range(t)
#         h.origin[1] + surface_vec[1] * width/2 + t * scale * surface_vec[1]
#     end
#
#     function y(t)
#         check_boundary_coord_range(t)
#         h.origin[2] + t * scale * surface_vec[2]
#     end
#
#     return x, y
# end
