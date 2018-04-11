"""
Shape with 4 right-angles and axis aligned sides, defined by width and height.
Note: the origin is in the center of the rectangle.
"""
struct Rectangle{T} <: Shape{2,T}
    origin::SVector{2,T}
    width::T
    height::T
end

# Alternative constructor, where bottomleft and topright corners are specified
function Rectangle(bottomleft::SVector{2,T}, topright::SVector{2,T}) where {T}
    origin = (bottomleft .+ topright)/2
    width = topright[1] - bottomleft[1]
    height = topright[2] - bottomleft[2]
    Rectangle{T}(origin, width, height)
end


name(shape::Rectangle) = "Rectangle"

volume(rectangle::Rectangle) = rectangle.width*rectangle.height

function congruent(r1::Rectangle{T}, r2::Rectangle{T}) where T
    r1.width  == r2.width  &&
    r1.height == r2.height
end

# Rectangle bounds itself
function bounding_rectangle(rect::Rectangle)
    rect
end

"Return SVector with the coordinates of the bottom left of a rectangle"
bottomleft(rect::Rectangle) = origin(rect) .- SVector(width/2, height/2)

"Return SVector with the coordinates of the top right of a rectangle"
topright(rect::Rectangle)   = origin(rect) .+ SVector(width/2, height/2)

# Create a box which bounds two shapes
function bounding_rectangle(shape1::Shape, shape2::Shape)
    box1 = bounding_rectangle(shape1)
    box2 = bounding_rectangle(shape2)

    min_bottomleft = min.(bottomleft(box1), bottomleft(box2))
    max_topright = max.(topright(box1), topright(box2))

    return Rectangle(min_bottomleft, max_topright)
end

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
