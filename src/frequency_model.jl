
type FrequencyModel{T <: AbstractFloat}
    shape::Shape
    particles::Vector{Particle{T}}
    response::Matrix{Complex{T}}
    hankel_order::Int
    k_arr::Vector{T}
    # Each column in this matrix is a vector at which we listen for the response
    listener_positions::Matrix{T}
    source_direction::Vector{T}
    seed::Vector{UInt32}
end

"""
Calcuate the volume fraction of this simulation from the volume of the particles
and the bounding shape.
"""
function volfrac(sim::FrequencyModel)
    volume(sim.particles)/volume(sim.shape)
end