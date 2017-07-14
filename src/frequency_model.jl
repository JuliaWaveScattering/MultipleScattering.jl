
type FrequencyModel{T <: AbstractFloat}
    shape::Shape
    ρ::T #the medium's density
    c::Complex{T} #the medium's phase velocity
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
function calculate_volfrac(model::FrequencyModel)
    volume(model.particles)/volume(model.shape)
end

"""
Find the mean radius of the particles in this model
"""
mean_radius(model::FrequencyModel) = mean_radius(model.particles)

"""
Find the standard deviation of the radii of the particles in this model
"""
std_radius(model::FrequencyModel) = std_radius(model.particles)
