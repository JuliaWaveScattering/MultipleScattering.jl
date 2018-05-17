"Shape where boundary is a fixed distance from the origin"
struct Circle{T} <: Shape{2,T}
    origin::SVector{2,T}
    radius::T
end

# Alternate constructors, where type is inferred naturally
Circle(origin::Tuple{T,T}, radius::T) where {T} = Circle{T}(origin, radius)
Circle(origin::Vector{T}, radius::T) where {T} = Circle{T}(origin, radius)

name(shape::Circle) = "Circle"

outer_radius(c::Circle) = c.radius
volume(shape::Circle) = π * shape.radius^2

function inside{T}(outer_circle::Circle{T}, inner_circle::Circle{T})
    norm(origin(outer_circle) - origin(inner_circle)) <= outer_circle.radius - inner_circle.radius
end

function inside(circle::Circle, x::AbstractVector)
    norm(origin(circle) .- x) <= circle.radius
end

function inside(rect::Rectangle{T}, circle::Circle{T}) where {T}
    all((origin(circle) .- circle.radius) .>= bottomleft(rect)) &&
    all((origin(circle) .+ circle.radius) .<= topright(rect))
end

function congruent(c1::Circle{T}, c2::Circle{T}) where T
    c1.radius == c2.radius
end

function bounding_rectangle{T}(circle::Circle{T})
    return Rectangle{T}(origin(circle) .- shape.radius, origin(circle) .+ shape.radius)
end

function boundary_functions{T}(circle::Circle{T})

    function x(t)
        check_boundary_coord_range(t)
        circle.radius * cos(2π * t) + origin(circle)[1]
    end

    function y(t)
        check_boundary_coord_range(t)
        circle.radius * sin(2π * t) + origin(circle)[2]
    end

    return x, y
end
