
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
    topright_particle(p) = p.x .+ p.r
    broadcastmax(x,y) = max.(x,y)
    topright = mapreduce(topright_particle, broadcastmax, particles)

    bottomleft_particle(p) = p.x .- p.r
    broadcastmin(x,y) = min.(x,y)
    bottomleft = mapreduce(bottomleft_particle, broadcastmin, particles)

    Rectangle(bottomleft, topright)
end

function inside{T}(shape::Rectangle{T}, particle::Particle{T})
    all(particle.x .- particle.r .> shape.bottomleft) &&
    all(particle.x .+ particle.r .< shape.topright)
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
    norm(shape.centre - particle.x) < shape.radius - particle.r
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
        cos(2π * t) + shape.centre[1]
    end
    
    function y(t) 
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end
        sin(2π * t) + shape.centre[2]
    end
    
    return x, y
end

# =============================== TimeOfFlight =================================
"""
We define a shape where all particles are less than time away from
listener_position, but also with a positive x coordinate (circle segment).
"""
type TimeOfFlight{T <: AbstractFloat} <: Shape
    listener_position::Vector{T}
    time::T
end

function inside{T}(shape::TimeOfFlight{T},particle::Particle{T})
    particle.x[1] > 0.0 &&
    norm(particle.x - shape.listener_position) < shape.time
end

function volume{T}(shape::TimeOfFlight{T})
    θ = 2acos(-shape.listener_position[1] / shape.time)
    return shape.time^2 * (θ - sin(θ)) / 2
end

function bounding_box{T}(shape::TimeOfFlight{T})
    box_height = 2sqrt(shape.time^2 - shape.listener_position[1]^2)
    box_width = max(shape.time + shape.listener_position[1], zero(T))
    return Rectangle([zero(T), -box_height / 2], [box_width, box_height / 2])
end

name(shape::TimeOfFlight) = "Time of flight"

function boundary_functions(shape::TimeOfFlight)
    function x(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end
        
        if t <= 1/2 
            θ = acos(-shape.listener_position[1] / shape.time)
            return shape.time*cos(θ*(4t-1)) + shape.listener_position[1]
        else
            return 0.0
        end
    end
    function y(t)
        if t<0 || t>1 error("Boundary coordinate must be between 0 and 1") end
        
        if t <= 1/2 
            θ = acos(-shape.listener_position[1] / shape.time)
            return shape.time*sin(θ*(4t-1)) + shape.listener_position[2]
        else
            return 2*(3/4-t)*2*sqrt(shape.time^2 - shape.listener_position[1]^2)
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
