abstract type Simulation{Dim,T} end

mutable struct FrequencySimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
end

TwoDimAcousticFrequencySimulation{T} = FrequencySimulation{2,Acoustic{2,T},T}

import Base.run

# Main run function, all other run functions use this
function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x::Vector{SVector{Dim,T}})::FrequencySimulationResult{Dim,P,T} where {Dim,P,T}
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

mutable struct TimeSimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
    impulse::Function
end
