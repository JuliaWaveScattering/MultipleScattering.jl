
"Abstract class for Results of Simulations"
abstract type SimulationResult{Dim,T} end

struct FrequencySimulationResult{Dim,FieldDim,P<:PhysicalProperties,T<:AbstractFloat} <: SimulationResult{Dim,T}
    field::Matrix{MVector{FieldDim,T}}
    x::Vector{MVector{Dim,T}}
    ω::RowVector{T}
end

struct TimeSimulationResult{Dim,FieldDim,P<:PhysicalProperties,T<:AbstractFloat} <: SimulationResult{Dim,T}
    field::Matrix{MVector{FieldDim,T}}
    x::Vector{MVector{Dim,T}}
    t::RowVector{T}
end


import Base.union
"""
Combine two FrequencyResults intelligently, allows user to optionally sort ω
and x
"""
function union(r1::FrequencySimulationResult, r2::FrequencySimulationResult; sort::Function = identity)
    if r1.x == r2.x
    elseif r1.k == r2.k
    end
end

function union(r1::SimulationResult,r2::SimulationResult)
    error("No implementation of union found for Simulation Results of type $(typeof(r1)) and $(typeof(r2))")
end
