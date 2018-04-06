abstract type Simulation{Dim,T} end

mutable struct FrequencySimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
end

TwoDimAcousticFrequencySimulation{T} = FrequencySimulation{2,Acoustic{2,T},T}

import Base.run

# Main run function, all other run functions use this
function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x::Vector{SVector{Dim,T}}) where {Dim,P,T}

    # Number of Hankel modes
    Nh = 5

    # Precompute T-matrices for these particles
    t_matrices = get_t_matrices(sim.medium, sim.particles, ω, Nh)

    # Compute scattering matrix
    S = scattering_matrix(sim.medium, sim.particles, t_matrices, ω, Nh)

    # Get forcing vector for this source
    f = forcing(sim.source, sim.particles, t_matrices, ω, Nh)

    # Find Hankel coefficients by solving scattering matrix for this forcing
    a = f\S

    # Evaluate Hankel series at the requested x positions

end

function run(sim::FrequencySimulation{Dim,P,T}, ω::Vector{T}, x::Vector{SVector{Dim,T}})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
    # Compute for each angular frequency, then join up all the results
    mapreduce(ω->run(sim,ω,x), union, ω)
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::Vector{T}, x::SVector{Dim,T})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
    run(sim,ω,[x])
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x::SVector{Dim,T})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
    run(sim,[ω],[x])
end

function forcing(source::Source{P,T}, particles::Vector{AbstractParticle{Dim,T}}, t_matrices::Vector{AbstractMatrix}, ω::T, Nh::Integer)::Vector{Complex{T}} where {Dim,P,T}
    mat = [source.coef(n,p.position,ω) for n in -Nh:Nh, p in particles]
    f = Vector{Complex{T}}(prod(size(mat)))
    H = 2Nh + 1
    for i in eachindex(particles)
        f[(i-1)*H+1:i*H] .= t_matrices[i]*mat[:,i]
    end
    return f
end

mutable struct TimeSimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
    impulse::Function
end
