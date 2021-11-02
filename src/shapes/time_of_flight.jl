"""
A shape which contains all particles,``x > x_0``, necessary to simulate a plane-wave scattering from an infinite medium, for a reciever at the focal point, for time ``t < D / c``, where ``c`` is the sound speed of the background medium, and ``D`` is some chosen focal distance. The plane-wave travels towards the positive direction of the ``x`` axis.

More precisely, if the focal point is at ``(x_f,y_f)`` then the interior of the shape
is defined as
``(y - y_f)^2 < (D + x_0)^2 - x_f^2 - 2(D + x_0 - x_f)x``  and ``x > min(x_0, x_f)``
where ``D`` is the focal distance.
"""
struct TimeOfFlightPlaneWaveToPoint{T <: AbstractFloat,Dim} <: Shape{T,Dim}
    focal_point::AbstractVector{T}
    focal_distance::T
    minimum_x::T
end

TimeOfFlightPlaneWaveToPoint(focal_point::AbstractVector{T},  focal_distance::T;  minimum_x::T = zero(T)) where T <:AbstractFloat = TimeOfFlightPlaneWaveToPoint{T,length(focal_point)}(focal_point, focal_distance, minimum_x)

name(shape::TimeOfFlightPlaneWaveToPoint) = "Time of flight from planar source to the focal point"

# NOTE: this needs to be redone
function volume(shape::TimeOfFlightPlaneWaveToPoint{T}) where T <: AbstractFloat

    D = shape.focal_distance
    xf = shape.focal_point[1]
    x_min = shape.minimum_x
    x_max = (x_min + xf + D) / T(2)

    a = D + x_min + xf
    b = D + x_min - xf
    return 2/3*(sqrt(b*(a - 2*x_min)^3) - sqrt(b*(a - 2*x_max)^3))
end

import Base.issubset
function issubset(sphere::Sphere{T,Dim}, shape::TimeOfFlightPlaneWaveToPoint{T,Dim}) where {T, Dim}
    return (origin(sphere)[1] > shape.minimum_x) && (norm(origin(sphere) - shape.focal_point) <= shape.focal_distance - origin(sphere)[1]  - T(2) * outer_radius(sphere))
end

function bounding_box(shape::TimeOfFlightPlaneWaveToPoint{T,Dim}) where {T,Dim}
    D = shape.focal_distance
    xf = shape.focal_point[1]
    y_f = shape.focal_point[2]
    x_min = shape.minimum_x
    x_max = (x_min + xf + D) / T(2)

    centre = zeros(Dim)
    centre[1] = (x_min + x_max) / T(2)
    centre[2] = y_f

    dimensions = ones(Dim) * T(2) * sqrt((D + x_min)^2 - xf^2 - 2*(D + x_min - xf)*x_min)
    dimensions[1] = x_max - x_min

    return Box(centre, dimensions)
end


function boundary_functions(shape::TimeOfFlightPlaneWaveToPoint{T,2}) where T

    D = shape.focal_distance
    xf = shape.focal_point[1]
    y_f = shape.focal_point[2]
    x_min = shape.minimum_x
    x_max = (x_min + xf + D) / T(2)

    function x(τ)
        check_boundary_coord_range(τ)

        if τ <= 1//3
            return (1 - 3*τ)*x_min + 3*τ*x_max
        elseif τ <= 2//3
            return (3*τ - 1)*x_min + (2 - 3*τ)*x_max
        else
            return x_min
        end
    end

    function y(τ)
        check_boundary_coord_range(τ)

        if τ <= 1//3
            return y_f + sqrt((D + x_min)^2 - xf^2 - 2*(D + x_min - xf)*((1 - 3*τ)*x_min + 3*τ*x_max))
        elseif τ <= 2//3
            return y_f - sqrt((D + x_min)^2 - xf^2 - 2*(D + x_min - xf)*((3*τ - 1)*x_min + 3*(2//3 - τ)*x_max))
        else
            return y_f + (6*τ - 5)*sqrt((D + x_min)^2 - xf^2 - 2*(D + x_min - xf)*x_min)
        end
    end

    return x, y
end
