
"Abstract class for Results of Simulations"
abstract type SimulationResult{Dim,T} end

struct FrequencySimulationResult{Dim,FieldDim,P<:PhysicalProperties,T<:AbstractFloat} <: SimulationResult{Dim,T}
    field::Matrix{SVector{FieldDim,Complex{T}}}
    x::Vector{SVector{Dim,T}}
    ω::RowVector{T}
end

struct TimeSimulationResult{Dim,FieldDim,P<:PhysicalProperties,T<:AbstractFloat} <: SimulationResult{Dim,T}
    field::Matrix{SVector{FieldDim,Complex{T}}}
    x::Vector{SVector{Dim,T}}
    t::RowVector{T}
end

"""
Get vector field from result as a matrix
"""
field(result::SimulationResult) = result.field

"""
Get vector field from result, with x-index i and ω-index j
"""
function field(result::FrequencySimulationResult, i::Integer, j::Integer)
    result.field[i,j]
end

"""
Get scalar field from result as a matrix
"""
function field(result::FrequencySimulationResult{Dim,1,P,T})::Matrix{Complex{T}} where {Dim, T, P<:PhysicalProperties{Dim,1,T}}
    map(x->x[1], result.field)
end

"""
Get scalar field from result, with x-index i and ω-index j
"""
function field(result::FrequencySimulationResult{Dim,1,P,T}, i::Integer, j::Integer)::Complex{T} where {Dim, T, P<:PhysicalProperties{Dim,1,T}}
    result.field[i,j][1]
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
