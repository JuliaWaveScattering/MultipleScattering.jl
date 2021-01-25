"""
    Rectangle([origin::AbstractVector{T}=zeros(),] width::T, Height::T)
    Rectangle(bottomleft::AbstractVector{T}, topright::AbstractVector{T})

2D [`Shape`](@ref) with axis aligned sides, defined by width, height and origin (at the center).
"""
struct Rectangle{T} <: Shape{T,2}
    origin::SVector{2,T}
    width::T
    height::T
end

# Alternative constructor, where bottomleft and topright corners are specified
function Rectangle(bottomleft::Union{AbstractVector{T},NTuple{2,T}}, topright::Union{AbstractVector{T},NTuple{2,T}}) where {T}
    origin = (bottomleft .+ topright)./2
    width = topright[1] - bottomleft[1]
    height = topright[2] - bottomleft[2]
    if sum(bottomleft .< topright) < 2
        error("bottomleft: $bottomleft should be below and left of topright: $topright.")
    end
    Rectangle{T}(origin, width, height)
end

# Alternate constructors, where type is inferred naturally
Rectangle(origin::Tuple{T,T}, width::T, height::T) where {T} = Rectangle{T}(origin, width, height)
Rectangle(origin::Vector{T}, width::T, height::T) where {T} = Rectangle{T}(origin, width, height)
# If no position is given, assume origin is at zero
Rectangle(width::T, height::T) where {T} = Rectangle{T}(SVector(zero(T),zero(T)), width, height)

name(shape::Rectangle) = "Rectangle"

outer_radius(r::Rectangle) = sqrt((r.width/2)^2 + (r.height/2)^2)
volume(rectangle::Rectangle) = rectangle.width*rectangle.height

import Base.issubset
function issubset(inner_rect::Rectangle{T}, outer_rect::Rectangle{T}) where T
    all(topright(inner_rect) .<= topright(outer_rect)) &&
    all(bottomleft(inner_rect) .>= bottomleft(outer_rect))
end

import Base.in
function in(x::AbstractVector, r::Rectangle)::Bool
    all(abs.(x .- r.origin) .<= SVector(r.width, r.height))
end

import Base.(==)
function ==(r1::Rectangle{T}, r2::Rectangle{T}) where T
    r1.origin == r2.origin &&
    r1.width  == r2.width  &&
    r1.height == r2.height
end

import Base.isequal
function isequal(r1::Rectangle{T}, r2::Rectangle{T}) where T
    isequal(r1.origin, r2.origin) &&
    isequal(r1.width , r2.width ) &&
    isequal(r1.height, r2.height)
end

function iscongruent(r1::Rectangle{T}, r2::Rectangle{T}) where T
    r1.width  == r2.width  &&
    r1.height == r2.height
end

function congruent(r::Rectangle{T}, x) where T
    Rectangle{T}(x, r.width, r.height)
end

# Rectangle bounds itself
bounding_box(rect::Rectangle) = rect

"Return SVector with the coordinates of the bottom left of a rectangle"
bottomleft(rect::Rectangle) = origin(rect) .- SVector(rect.width/2, rect.height/2)

"Return SVector with the coordinates of the top right of a rectangle"
topright(rect::Rectangle)   = origin(rect) .+ SVector(rect.width/2, rect.height/2)


function boundary_functions(rect::Rectangle)

    function x(t)
        check_boundary_coord_range(t)
        if     t <= 1//4 x = bottomleft(rect)[1] + 4t*rect.width
        elseif t <= 2//4 x = topright(rect)[1]
        elseif t <= 3//4 x = topright(rect)[1] - 4*(t-2//4)*rect.width
        else             x = bottomleft(rect)[1]
        end
    end

    function y(t)
        check_boundary_coord_range(t)
        if     t <= 1//4 x = bottomleft(rect)[2]
        elseif t <= 2//4 x = bottomleft(rect)[2] + 4*(t-1//4)*rect.height
        elseif t <= 3//4 x = topright(rect)[2]
        else             x = topright(rect)[2] - 4*(t-3//4)*rect.height
        end
    end

    return x, y
end
