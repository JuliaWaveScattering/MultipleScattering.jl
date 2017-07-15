
# ================================ Rectangle ==================================
# We define a rectangle using the two corners where topright .> bottomleft
type Rectangle{T <: AbstractFloat} <: Shape
    bottomleft::Vector{T}
    topright::Vector{T}
end

"Generates a rectangle for Specified volfrac, radius and num_particles"
function Rectangle{T}(volfrac::Number, radius::T, num_particles::Int)
    w = sqrt(num_particles*radius^2*pi/volfrac)
    Rectangle([0.0,-w/2], [w,w/2])
end

"Generates a rectangle which contains all the particles"
function Rectangle{T}(particles::Vector{Particle{T}})
    temp_max(M) = maximum(M, 2)
    topright_particle(p) = p.x .+ p.r
    topright = mapreduce(topright_particle, temp_max, particles)

    bottomleft_particle(p) = p.x .- p.r
    temp_min(M) = minimum(M, 2)
    bottomleft = mapreduce(bottomleft_particle, temp_min, particles)

    Rectangle(bottomleft, topright)
end

function inside{T}(shape::Rectangle{T}, particle::Particle{T})
    all(particle.x .- particle.r .> shape.bottomleft) && all(particle.x .+ particle.r .< shape.topright)
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
        if 0 <= t <= 1/4 x = shape.bottomleft[1]
        elseif t <= 2/4 x = shape.bottomleft[1]
        elseif t <= 3/4 x = shape.bottomleft[1]
        elseif t <= 4/4 x = shape.bottomleft[1]
        else error("Boundary coordinate t should only be between 0 and 1")
        end
    end
    function y(t)
        if 0 <= t <= 1/4 x = shape.bottomleft[1]
        elseif t <= 2/4 x = shape.bottomleft[1]
        elseif t <= 3/4 x = shape.bottomleft[1]
        elseif t <= 4/4 x = shape.bottomleft[1]
        else error("Boundary coordinate t should only be between 0 and 1")
        end
    end
    return (x, y)
end


# ================================= Circle ===================================
# Circle defined by a centre and a radius
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
    x(t) = shape.cos(2π * t) + shape.centre[1]
    y(t) = shape.sin(2π * t) + shape.centre[2]
    return (x, y)
end

# =============================== TimeOfFlight =================================
# We define a shape where all particles are less than time away from
# listener_position, but also with a positive x coordinate (circle segment).
type TimeOfFlight{T <: AbstractFloat} <: Shape
    listener_position::Vector{T}
    time::T
end

function inside{T}(shape::TimeOfFlight{T},particle::Particle{T})
    particle.x[1] > 0.0 && norm(particle.x - shape.listener_position) < shape.time
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

name(shape::TimeOfFlight) = "TimeOfFlight"

function boundary_functions(shape::TimeOfFlight)
    function x(t)
        if 0 <= t <= 1/2 x = 0
        elseif t <= 2/2 x = bottomleft[1]
        else error("Boundary coordinate t should only be between 0 and 1")
        end
    end
    function y(t)
        if 0 <= t <= 1/2 x = bottomleft[1]
        elseif t <= 2/2 x = bottomleft[1]
        else error("Boundary coordinate t should only be between 0 and 1")
        end
    end
    return (x, y)
end
