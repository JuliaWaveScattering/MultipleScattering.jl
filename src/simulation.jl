abstract type Simulation{T,Dim} end

mutable struct FrequencySimulation{P,T,Dim} <: Simulation{T,Dim} where T <: AbstractFloat, P <: PhysicalProperties{Dim,FieldDim,T}, Dim::Int
    medium::P
    particles::Vector{Particle{P,T,Dim}}
    source::Source{P,T}
end

# Type aliases for convenience
TwoDimAcousticFrequencySimulation{T} = FrequencySimulation{Acoustics{T},T,2}


function run(sim::FrequencySimulation{P,T,Dim}, ω::T, x::SVector{Dim,T})
end

function run(sim::FrequencySimulation{P,T,Dim}, ω::T, x::Vector{SVector{Dim,T}})
end

function run(sim::FrequencySimulation{P,T,Dim}, ω::Vector{T}, x::SVector{Dim,T})
end

function run(sim::FrequencySimulation{P,T,Dim}, ω::Vector{T}, x::Vector{SVector{Dim,T}})
end



mutable struct TimeSimulation{P,T,Dim} <: Simulation{T,Dim} where T <: AbstractFloat, P <: PhysicalProperties{Dim,FieldDim,T}, Dim::Int
    medium::P
    particles::Vector{Particle{P,T,Dim}}
    source::Source{P,T}
    impulse::Impulse{P,T}
end
