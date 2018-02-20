
# ================================ Rectangle ==================================
"We define a rectangle using the two corners where topright .> bottomleft"
type Rectangle{T <: AbstractFloat} <: Shape
    bottomleft::Vector{T}
    topright::Vector{T}
end

"Generates a rectangle for Specified volfrac, radius and num_particles"
function Rectangle{T}(volfrac::Number, radius::T, num_particles::Number)
    w = sqrt(num_particles*radius^2*pi/volfrac)
    Rectangle([zero(T),-w/2], [w,w/2])
end

"Generates a rectangle which contains all the particles"
function Rectangle{T}(particles::Vector{Particle{T}})

    if isempty(particles) return Rectangle([zero(T),zero(T)], [zero(T),zero(T)]) end
    topright_particle(p) = p.x .+ p.r
    broadcastmax(x,y) = max.(x,y)
    topright = mapreduce(topright_particle, broadcastmax, particles)

    bottomleft_particle(p) = p.x .- p.r
    broadcastmin(x,y) = min.(x,y)
    bottomleft = mapreduce(bottomleft_particle, broadcastmin, particles)

    Rectangle(bottomleft, topright)
end

function inside{T}(shape::Rectangle{T}, particle::Particle{T})
    all(particle.x .- particle.r .>= shape.bottomleft) &&
    all(particle.x .+ particle.r .<= shape.topright)
end

function volume{T}(shape::Rectangle{T})
    return prod(shape.topright .- shape.bottomleft)
end

# Rectangles already are bounding boxes
function bounding_box{T}(shape::Rectangle{T})
    shape
end

name(shape::Rectangle) = "Rectangle"

function boundary_functions(shape::Rectangle)
    function x(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end

        width = shape.topright[1] - shape.bottomleft[1]
        if     t <= 1/4 x = shape.bottomleft[1] + 4t*width
        elseif t <= 2/4 x = shape.topright[1]
        elseif t <= 3/4 x = shape.topright[1] - 4*(t-2/4)*width
        else            x = shape.bottomleft[1]
        end
    end
    function y(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end

        height = shape.topright[2] - shape.bottomleft[2]
        if     t <= 1/4 x = shape.bottomleft[2]
        elseif t <= 2/4 x = shape.bottomleft[2] + 4*(t-1/4)*height
        elseif t <= 3/4 x = shape.topright[2]
        else            x = shape.topright[2] - 4*(t-3/4)*height
        end
    end
    return x, y
end


# ================================= Circle ===================================
"Circle defined by a centre and a radius"
type Circle{T <: AbstractFloat} <: Shape
    radius::T
    centre::Vector{T}
end

function inside{T}(shape::Circle{T}, particle::Particle{T})
    norm(shape.centre - particle.x) <= shape.radius - particle.r
end

function volume{T}(shape::Circle{T})
    return π * shape.radius^2
end

function bounding_box{T}(shape::Circle{T})
    return Rectangle(shape.centre .- shape.radius, shape.centre .+ shape.radius)
end

name(shape::Circle) = "Circle"

function boundary_functions{T}(shape::Circle{T})
    function x(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end
        shape.radius * cos(2π * t) + shape.centre[1]
    end

    function y(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end
        shape.radius * sin(2π * t) + shape.centre[2]
    end

    return x, y
end

# =============================== TimeOfFlight =================================
"""
A shape where anything inside could cause a disturbance at the listener position
from a planar wavefront parallel to the y axis starting at the listener. Also
everything inside has a positive `x` coordinate.

More precisely, if the listener is at (l_x,l_y) then the interior of the shape
is defined as
x-l_x+sqrt((x-l_x)^2+(y-l_y)^2)<time and x>0
"""
type TimeOfFlight{T <: AbstractFloat} <: Shape
    listener_position::Vector{T}
    time::T
end

function inside{T}(shape::TimeOfFlight{T},particle::Particle{T})
    l_to_p = particle.x - shape.listener_position
    return (particle.x[1] > 0.0) && (l_to_p[1] + norm(l_to_p) <= shape.time)
end

# Let T = shape.time, xr = shape.listener_position[1], yr = shape.listener_position[2]
# and assume all particles are placed in x>0.
# Let (x,y) be a point on the curved part of the shape, then
# x - xr + sqrt((x - xr)^2 + (y - yr)^2) == T => x == T/2 + xr - (y - yr)^2/(2T)
# then the area = 2*Integrate[ T/2 + xr - (y - yr)^2/(2T), {y,yr, yr + sqrt(T^2 + 2xr*T)}]
function volume{T}(shape::TimeOfFlight{T})
    l_x = shape.listener_position[1]
    return 2/(3*shape.time)*(shape.time^2 + 2*l_x*shape.time)^(3//2)
end

function bounding_box{T}(shape::TimeOfFlight{T})
    t = shape.time
    l = shape.listener_position
    x_max = max(t/2 + l[1], zero(T))
    return Rectangle([ zero(T), l[2] - sqrt(t^2 + 2t*l[1]) ],
                     [ x_max,   l[2] + sqrt(t^2 + 2t*l[1]) ])
end

name(shape::TimeOfFlight) = "Time of flight"

function boundary_functions{T}(shape::TimeOfFlight{T})
    function x(τ)
        if τ<0 || τ>1 error("Boundary coordinate must be between 0 and 1") end
        if τ <= 1//2
            return zero(T)
        else
            t = shape.time
            l = shape.listener_position
            y_sq = (4*(3//4-τ))^2*(t^2 + 2t*l[1])
            return l[1] + (t^2-y_sq)/(2t)
        end
    end
    function y(τ)
        if τ<0 || τ>1 error("Boundary coordinate must be between 0 and 1") end
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

# ========================== TimeOfFlightFromPoint ============================
"""
A shape where anything inside could cause a disturbance at the listener position
from a point source wavefront starting at the listener. Also everything inside
has a positive `x` coordinate. It is equivalent to a segment of a circle.

More precisely, if the listener is at (l_x,l_y) then the interior of the shape
is defined as
sqrt((x-l_x)^2+(y-l_y)^2)<time and x>0
"""
type TimeOfFlightFromPoint{T <: AbstractFloat} <: Shape
    listener_position::Vector{T}
    time::T
end

function inside{T}(shape::TimeOfFlightFromPoint{T},particle::Particle{T})
    particle.x[1] > 0 &&
    norm(particle.x - shape.listener_position) < shape.time
end

function volume{T}(shape::TimeOfFlightFromPoint{T})
    θ = 2acos(-shape.listener_position[1] / shape.time)
    return shape.time^2 * (θ - sin(θ)) / 2
end

function bounding_box{T}(shape::TimeOfFlightFromPoint{T})
    box_height = 2sqrt(shape.time^2 - shape.listener_position[1]^2)
    box_width = max(shape.time + shape.listener_position[1], zero(T))
    return Rectangle([zero(T), -box_height / 2], [box_width, box_height / 2])
end

name(shape::TimeOfFlightFromPoint) = "Time of flight from point"

function boundary_functions(shape::TimeOfFlightFromPoint)
    function x(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end

        if t <= 1//2
            θ = acos(-shape.listener_position[1] / shape.time)
            return shape.time*cos(θ*(4t-1)) + shape.listener_position[1]
        else
            return 0.0
        end
    end
    function y(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end

        if t <= 1//2
            θ = acos(-shape.listener_position[1] / shape.time)
            return shape.time*sin(θ*(4t-1)) + shape.listener_position[2]
        else
            return 2*(3//4-t)*2*sqrt(shape.time^2 - shape.listener_position[1]^2)
        end
    end
    return x, y
end

# =============================== Utility =================================
"Create a box which bounds two shapes"
function bounding_box(shape1::Shape,shape2::Shape)
    box1 = bounding_box(shape1)
    box2 = bounding_box(shape2)

    bottomleft = min.(box1.bottomleft,box2.bottomleft)
    topright = max.(box1.topright,box2.topright)

    return Rectangle(bottomleft,topright)
end

"Generates a rectangle which contains all the particles"
function bounding_box{T}(particles::Vector{Particle{T}})
    Rectangle(particles)
end
