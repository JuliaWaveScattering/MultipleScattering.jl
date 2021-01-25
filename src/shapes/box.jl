"""
    Box([origin::AbstractVector{T}=zeros(),] dimensions::AbstractVector{T})

A [`Box`](@ref) for 2D and 3D with axis aligned sides, defined by dimensions and origin (at the center).
"""
struct Box{T,Dim} <: Shape{T,Dim}
    origin::SVector{Dim,T}
    dimensions::SVector{Dim,T}
end

# Alternative constructor, where bottomleft and topright corners are specified
function Box(origin::AbstractVector{T}, dimensions::AbstractVector{T}) where T
    Dim = length(origin)
    Box{T,Dim}(origin, dimensions)
end

function Box(origin::NTuple{Dim,T}, dimensions::NTuple{Dim,T}) where {T,Dim}
    Box{T,Dim}(origin, dimensions)
end

function Box(corners::Vector{S}) where S<:AbstractVector
    centre = mean(corners)
    cs = hcat([abs.(c - centre) for c in corners]...)
    dimensions = 2 .* abs.(corners[1] - centre)

    if dimensions != 2 .* abs.(corners[2] - centre)
        @error "expected $corners to be a vector of corners of a $(length(corners[1])) dimensional box, with sides aligned with the coordinate axis."
    end

    Box(centre,dimensions)
end

Box(dimensions::AbstractVector{T}) where T = Box(zeros(T,length(dimensions)),dimensions)

name(shape::Box) = "Box"
name(shape::Box{T,2}) where T = "Rectangle"

outer_radius(r::Box) = sqrt(sum(r.dimensions .^ 2))
volume(box::Box) = prod(box.dimensions)

import Base.issubset
function issubset(inner::Box, outer::Box)
    corners = box_corners(inner)
    return all(c ∈ outer for c in corners)
end

import Base.in
function in(x::AbstractVector, b::Box{T})::Bool where T
    all(abs.(x .- b.origin) .<= b.dimensions ./ T(2))
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
function box_corners(b::Box{T,Dim}) where {T,Dim}
    d1 = origin(b) - b.dimensions ./ T(2)
    d2 = origin(b) + b.dimensions ./ T(2)
    ds =[d1,d2]

    corners = map(0:(2^Dim-1)) do i
        switches = 1 .+ [parse(Int,c)  for c in bitstring(i)[end-Dim+1:end]]
        [ds[switches[j]][j] for j in eachindex(switches)]
    end

    return corners
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
        if     t <= 1//4 x = bottomleft(rect)[2]
        elseif t <= 2//4 x = bottomleft(rect)[2] + 4*(t-1//4) * h
        elseif t <= 3//4 x = topright(rect)[2]
        else             x = topright(rect)[2] - 4*(t-3//4) * h
        end
    end

    return x, y
end
