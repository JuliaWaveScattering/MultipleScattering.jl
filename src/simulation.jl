abstract type Simulation{Dim,T} end

mutable struct FrequencySimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
end

# Type aliases for convenience
if !isdefined(:TwoDimAcousticFrequencySimulation)
    TwoDimAcousticFrequencySimulation{T} = FrequencySimulation{2,Acoustics{2,T},T}
end

import Base.run

function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x::MVector{Dim,T}) where {Dim,P,T}
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::T, x::Vector{MVector{Dim,T}}) where {Dim,P,T}
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::Vector{T}, x::MVector{Dim,T}) where {Dim,P,T}
end

function run(sim::FrequencySimulation{Dim,P,T}, ω::Vector{T}, x::Vector{MVector{Dim,T}}) where {Dim,P,T}
end

mutable struct TimeSimulation{Dim,P<:PhysicalProperties,T<:AbstractFloat} <: Simulation{Dim,T}
    medium::P
    particles::Vector{AbstractParticle{Dim,T}}
    source::Source{P,T}
    impulse::Function
end
