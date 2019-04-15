"""
A shape where anything inside could cause a disturbance at the listener position
from a planar wavefront parallel to the y axis starting at the listener. Also
everything inside has a positive `x` coordinate.

More precisely, if the listener is at (l_x,l_y) then the interior of the shape
is defined as
x-l_x+sqrt((x-l_x)^2+(y-l_y)^2)<time and x>0
"""
struct TimeOfFlight{T <: AbstractFloat} <: Shape{T,2}
    listener_position::Vector{T}
    time::T
end

TimeOfFlight(pos::AbstractVector{T}, time::T) where T <:AbstractFloat = TimeOfFlight(Vector{T}(pos), time)

name(shape::TimeOfFlight) = "Time of flight from planar source"

# Let T = shape.time, xr = shape.listener_position[1], yr = shape.listener_position[2]
# and assume all particles are placed in x>0.
# Let (x,y) be a point on the curved part of the shape, then
# x - xr + sqrt((x - xr)^2 + (y - yr)^2) == T => x == T/2 + xr - (y - yr)^2/(2T)
# then the area = 2*Integrate[ T/2 + xr - (y - yr)^2/(2T), {y,yr, yr + sqrt(T^2 + 2xr*T)}]
function volume(shape::TimeOfFlight{T}) where T <: AbstractFloat
    l_x = shape.listener_position[1]
    return 2/(3*shape.time)*(shape.time^2 + 2*l_x*shape.time)^(3//2)
end

import Base.issubset
function issubset(circle::Circle, shape::TimeOfFlight)
    l_to_p = origin(circle) - shape.listener_position
    return (origin(circle)[1] > 0) && (l_to_p[1] + norm(l_to_p) <= (shape.time - 2circle.radius))
end

function bounding_rectangle(shape::TimeOfFlight{T}) where T <: AbstractFloat
    t = shape.time
    l = shape.listener_position
    x_max = max(t/2 + l[1], zero(T))
    return Rectangle(SVector( zero(T), l[2] - sqrt(t^2 + 2t*l[1])),
                     SVector( x_max,   l[2] + sqrt(t^2 + 2t*l[1])) )
end


function boundary_functions(shape::TimeOfFlight)

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
