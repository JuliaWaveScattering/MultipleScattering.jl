"""
A shape which contains all particles, with ``x > x_0``, necessary to simulate a plane-wave scattering from an infinite medium, for a reciever at the focal point, for time ``t < D / c``, where ``c`` is the sound speed of the background medium, and ``D`` is some chosen focal distance.

More precisely, if the focal point is at ``(x_f,y_f)`` then the interior of the shape
is defined as
``y^2 < D^2 + 2(x_f - D)x - x_f^2``  and ``x > x_0``
where ``D`` is the focal distance.
"""
struct TimeOfFlightPlaneWaveToPoint{T <: AbstractFloat,Dim} <: Shape{T,Dim}
    focal_point::Vector{T}
    focal_distance::T
    minimum_x::T
end

TimeOfFlightPlaneWaveToPoint(pos::AbstractVector{T}, time::T) where T <:AbstractFloat = TimeOfFlightPlaneWaveToPoint(Vector{T}(pos), time)

name(shape::TimeOfFlightPlaneWaveToPoint) = "Time of flight from planar source to the focal point"

# NOTE: this needs to be redone
function volume(shape::TimeOfFlightPlaneWaveToPoint{T}) where T <: AbstractFloat
    l_x = shape.listener_position[1]
    return 2/(3*shape.time)*(shape.time^2 + 2*l_x*shape.time)^(3//2)
end

import Base.issubset
function issubset(sphere::Sphere{T,Dim}, shape::TimeOfFlightPlaneWaveToPoint{T,Dim}) where {T, Dim}
    return (origin(sphere)[1] > shape.minimum_x) && (norm(origin(sphere) - shape.focal_point) <= shape.focal_distance - origin(sphere)[1]  - T(2) * outer_radius(sphere))
end

function bounding_box(shape::TimeOfFlightPlaneWaveToPoint{T,Dim}) where {T,Dim}
    D = shape.focal_distance
    xf = shape.focal_point
    x_min = shape.minimum_x
    x_max = (x_min + D) / T(2)

    centre = zeros(Dim)
    centre[1] = (x_min + x_max) / T(2)

    dimensions = ones(Dim) * T(2) * sqrt(D^2 + 2 * (xf - D) * x0 - xf^2 )
    dimensions[1] = x_max - x_min

    return Box(centre, dimensions)
end


function boundary_functions(shape::TimeOfFlightPlaneWaveToPoint{T,2}) where T

    function x(τ)
        check_boundary_coord_range(τ)
        if τ <= 1//2
            return zero(τ)
        else
            t = shape.time
            l = shape.listener_position
            y_sq = (4*(3//4-τ))^2*(t^2 + 2t*l[1])
            return l[1] + (t^2-y_sq)/(2t)
        end
    end

    function y(τ)
        check_boundary_coord_range(τ)
        t = shape.time
        l = shape.listener_position
        if τ <= 1//2
            return l[2] + 4*(τ-1//4)*sqrt(t^2 + 2t*l[1])
        else
            return l[2] + 4*(3//4-τ)*sqrt(t^2 + 2t*l[1])
        end
    end

    return x, y
end
