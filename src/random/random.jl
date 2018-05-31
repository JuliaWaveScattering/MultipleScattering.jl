
"Distribution of multiple scattering simulations"
abstract type SimulationDistribution{T,Dim} end

"""
Sample a simulation distribution returning a simultion. If no seed is provided,
a seed is generated.
"""
function sample(dist::SimulationDistribution, seed)
    error("Function sample is not implemented for this simulation distribution")
end

function sample(dist::SimulationDistribution)
    seed = generate_seed()
    sample(dist, seed)
end


"""
A distribution of frequency simulations where particles are all identical but
with random positions, contained within some specified shape with a specified
density. We can sample this to create a runnable FrequencySimulation.
"""
mutable struct FrequencySimulationDistribution{Dim,P,S,PP,PS,T} <: SimulationDistribution{T,Dim}
    medium::P
    shape::S
    volfrac::T
    source::Source{P,T}
    particle_medium::PP
    particle_shape::PS
end

function sample(dist::FrequencySimulationDistribution{Dim,P,S,PP,PS,T}, seed)::FrequencySimulation{T,Dim,P} where {Dim,P,S,PP,PS,T}

    medium = dist.medium
    particles = generate_particles(dist.shape, dist.volfrac, dist.particle_shape, dist.particle_medium, seed)
    source = dist.source

    return FrequencySimulation{T,Dim,P}(medium,particles,source)
end
