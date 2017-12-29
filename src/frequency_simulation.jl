
type FrequencySimulation{T <: AbstractFloat}
    shape::Shape
    ρ::T #the medium's density
    c::Complex{T} #the medium's phase velocity
    volfrac::T # the volume fraction of the particles
    particles::Vector{Particle{T}}
    response::Matrix{Complex{T}}
    hankel_order::Int
    k_arr::Vector{T}
    # Each column in this matrix is a vector at which we listen for the response
    listener_positions::Matrix{T}
    source_position::Vector{T}
    source_direction::Vector{T}
    seed::Vector{UInt32}
end

"""
Calcuate the volume fraction of this simulation from the volume of the particles
and the bounding shape.
"""
function calculate_volfrac(simulation::FrequencySimulation)
    volume(simulation.particles)/volume(simulation.shape)
end

"""
Find the mean radius of the particles in this simulation
"""
mean_radius(simulation::FrequencySimulation) = mean_radius(simulation.particles)

"""
Find the standard deviation of the radii of the particles in this simulation
"""
std_radius(simulation::FrequencySimulation) = std_radius(simulation.particles)
