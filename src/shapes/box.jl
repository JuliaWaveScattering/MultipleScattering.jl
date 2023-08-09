"""
    Box([origin::AbstractVector{T}=zeros(),] dimensions::AbstractVector{T})

A [`Box`](@ref) for 2D and 3D with axis aligned sides, defined by dimensions and origin (at the center).
"""
struct Box{T,Dim} <: Shape{Dim}
    origin::SVector{Dim,T}
    dimensions::SVector{Dim,T}
end

function Box(origin::AbstractVector{T}, dimensions::AbstractVector{T}) where T
    Dim = length(origin)
    Box{T,Dim}(origin, dimensions)
end

function Box(origin::NTuple{Dim,T}, dimensions::NTuple{Dim,T}) where {Dim,T}
    Box{T,Dim}(origin, dimensions)
end

Box(dimensions::AbstractVector{T}) where T = Box(zeros(T,length(dimensions)),dimensions)

"""
    Box(points::Vector{v} where v <: AbstractVector)

A [`Box`](@ref) for any dimension with axis aligned sides, that is a minimal cover for the points.
"""
function Box(points::Vector{v} where v <: AbstractVector) 
    ind = CartesianIndices(points[1])
    xs = [p[i] for p in points, i in ind]

    xmin = minimum(xs; dims = 1)
    xmax = maximum(xs; dims = 1)

    c = (xmin + xmax)[:] ./ 2.0;
    dims = (xmax - xmin)[:];

    return Box(c,dims)
end

Rectangle(bottomleft::Union{AbstractVector{T},NTuple{2,T}}, topright::Union{AbstractVector{T},NTuple{2,T}}) where T = Box([bottomleft,topright])

function Shape(b::Box{T,Dim};
        addtodimensions = zeros(T,Dim),
        vector_translation::AbstractVector{T} = zeros(T,Dim)
    ) where {T,Dim}
    Box(b.origin + vector_translation, b.dimensions .+ addtodimensions)
end

name(shape::Box) = "Box"
name(shape::Box{T,2}) where T = "Rectangle"

outer_radius(r::Box) = sqrt(sum(r.dimensions .^ 2))
volume(box::Box) = prod(box.dimensions)

import Base.issubset
function issubset(inner::Box, outer::Box)
    boxcorners = corners(inner)
    return all(c ∈ outer for c in boxcorners)
end

import Base.in
function in(x::AbstractVector, b::Box)::Bool
    all(abs.(x .- b.origin) .<= b.dimensions ./ 2)
end

import Base.(==)
function ==(r1::Box, r2::Box)
    r1.origin == r2.origin &&
    r1.dimensions  == r2.dimensions
end

import Base.isequal
function isequal(r1::Box, r2::Box)
    isequal(r1.origin, r2.origin) &&
    isequal(r1.dimensions , r2.dimensions )
end

function iscongruent(r1::Box, r2::Box)
    r1.dimensions  == r2.dimensions
end

function congruent(r::Box, x)
    Box(x, r.dimensions)
end

# Box bounds itself
bounding_box(b::Box) = b

"Returns a vector of SVector with the coordinates of the corners of the box"
function corners(b::Box{T,Dim}) where {T,Dim}
    d1 = origin(b) - b.dimensions ./ T(2)
    d2 = origin(b) + b.dimensions ./ T(2)
    ds =[d1,d2]

    boxcorners = map(0:(2^Dim-1)) do i
        switches = 1 .+ [parse(Int,c)  for c in bitstring(i)[end-Dim+1:end]]
        [ds[switches[j]][j] for j in eachindex(switches)]
    end

    return boxcorners
end

"Return SVector with the coordinates of the bottom left of a rectangle"
bottomleft(rect::Box{T,2}) where T = origin(rect) .- SVector(rect.dimensions) ./ 2

"Return SVector with the coordinates of the top right of a rectangle"
topright(rect::Box{T,2}) where T = origin(rect) .+ SVector(rect.dimensions) ./ 2


function boundary_functions(rect::Box{T,2}) where T
    w, h = rect.dimensions
    function x(t)
        check_boundary_coord_range(t)
        if     t <= 1//4 x = bottomleft(rect)[1] + 4t * w
        elseif t <= 2//4 x = topright(rect)[1]
        elseif t <= 3//4 x = topright(rect)[1] - 4*(t-2//4) * w
        else             x = bottomleft(rect)[1]
        end
    end

    function y(t)
        check_boundary_coord_range(t)
        if     t <= 1//4 y = bottomleft(rect)[2]
        elseif t <= 2//4 y = bottomleft(rect)[2] + 4*(t-1//4) * h
        elseif t <= 3//4 y = topright(rect)[2]
        else             y = topright(rect)[2] - 4*(t-3//4) * h
        end
    end

    return x, y
end
