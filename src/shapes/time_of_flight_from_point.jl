"""
A shape where anything inside could cause a disturbance at the listener position
from a point source wavefront starting at the listener. Also everything inside
has a positive `x` coordinate. It is equivalent to a segment of a circle.

More precisely, if the listener is at (l_x,l_y) then the interior of the shape
is defined as
sqrt((x-l_x)^2+(y-l_y)^2)<time and x>0
"""
struct TimeOfFlightPointWaveToPoint{T <: AbstractFloat} <: Shape{2}
    listener_position::Vector{T}
    time::T
end

name(shape::TimeOfFlightPointWaveToPoint) = "Time of flight from point source"

function volume(shape::TimeOfFlightPointWaveToPoint)
    θ = 2acos(-shape.listener_position[1] / shape.time)
    return shape.time^2 * (θ - sin(θ))/2
end

import Base.issubset
function issubset(circle::Sphere{2}, shape::TimeOfFlightPointWaveToPoint)
    (origin(circle)[1] - circle.radius) > 0 &&
    norm(origin(circle) - shape.listener_position) < (shape.time - 2circle.radius)
end

function bounding_box(shape::TimeOfFlightPointWaveToPoint)
    box_height = 2sqrt(shape.time^2 - shape.listener_position[1]^2)
    box_width = max(shape.time + shape.listener_position[1], 0)
    return Box([SVector(zero(box_width), -box_height/2), SVector(box_width, box_height/2)])
end


function boundary_functions(shape::TimeOfFlightPointWaveToPoint)

    function x(t)
        check_boundary_coord_range(t)
        if t <= 1//2
            θ = acos(-shape.listener_position[1] / shape.time)
            return shape.time*cos(θ*(4t-1)) + shape.listener_position[1]
        else
            return zero(t)
        end
    end

    function y(t)
        check_boundary_coord_range(t)
        if t <= 1//2
            θ = acos(-shape.listener_position[1] / shape.time)
            return shape.time*sin(θ*(4t-1)) + shape.listener_position[2]
        else
            return 2*(3//4-t)*2*sqrt(shape.time^2 - shape.listener_position[1]^2)
        end
    end

    return x, y
end
