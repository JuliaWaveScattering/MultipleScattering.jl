
mutable struct RandomFrequencySimulation{P,T,Dim} <: Simulation{T,Dim} where T <: AbstractFloat, P <: PhysicalProperties{T,Dim,FieldDim}, Dim::Int
    medium::P
    source::Source{P,T}
    shape::Shape
    volfrac::T
    radius::T
end

function run(sim::RandomFrequencySimulation, ω, x)
    seed = generate_seed()
    run(sim, ω, x, seed)
end

function run(sim::RandomFrequencySimulation{P,T,Dim}, ω::T, x::SVector{Dim,T}, seed)
end

function run(sim::RandomFrequencySimulation{P,T,Dim}, ω::T, x::Vector{SVector{Dim,T}}, seed)
end

function run(sim::RandomFrequencySimulation{P,T,Dim}, ω::Vector{T}, x::SVector{Dim,T}, seed)
end

function run(sim::RandomFrequencySimulation{P,T,Dim}, ω::Vector{T}, x::Vector{SVector{Dim,T}}, seed)
end

struct RandomFrequencySimulationResult{P,T,Dim} <: Result where P <: PhysicalProperties{T,Dim,FieldDim}, T <: AbstractFloat, Dim::Int
    field::Matrix{MVector{T,FieldDim}}
    x::Vector{MVector{T,FieldDim}}
    ω::RowVector{T}
    seed::Vector{UInt32}
end
