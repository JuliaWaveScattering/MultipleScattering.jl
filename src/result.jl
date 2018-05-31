
"Abstract class for Results of Simulations"
abstract type SimulationResult{T,Dim,FieldDim} end

struct FrequencySimulationResult{T<:AbstractFloat,Dim,FieldDim} <: SimulationResult{T,Dim,FieldDim}
    field::Matrix{SVector{FieldDim,Complex{T}}}
    x::Vector{SVector{Dim,T}}
    ω::RowVector{T}
end

function FrequencySimulationResult(field::Matrix{SVector{FieldDim,Complex{T}}}, x::AbstractVector{SVector{Dim,T}}, ω::AbstractVector{T}) where {Dim,T,FieldDim}
    FrequencySimulation{T,Dim,FieldDim}(field, Vector(x), RowVector(ω))
end

struct TimeSimulationResult{T<:AbstractFloat,Dim,FieldDim} <: SimulationResult{T,Dim,FieldDim}
    field::Matrix{SVector{FieldDim,T}}
    x::Vector{SVector{Dim,T}}
    t::RowVector{T}
end

function TimeSimulationResult(time_field::Union{Matrix{T},Matrix{AbstractVector{T}}}, x::AbstractVector{SVector{Dim,T}}, t::AbstractVector{T}) where {Dim,T}
    time_field = [SVector(d...) for d in time_field]
    FieldDim = size(time_field[1],1)
    TimeSimulationResult{T,Dim,FieldDim}(time_field, Vector(x), RowVector(t))
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
function field(result::FrequencySimulationResult{T,Dim,1})::Matrix{Complex{T}} where {Dim, T}
    map(x->x[1], result.field)
end

"""
Get scalar field from result, with x-index i and ω-index j
"""
function field(result::FrequencySimulationResult{T,Dim,1}, i::Integer, j::Integer)::Complex{T} where {Dim, T}
    result.field[i,j][1]
end

import Base.union
"""
Combine two FrequencyResults intelligently, allows user to optionally sort ω
and x
"""
function union(r1::FrequencySimulationResult, r2::FrequencySimulationResult; sort::Function = identity)
    if r1.x == r2.x
    elseif r1.ω == r2.ω
    end
end

function union(r1::SimulationResult,r2::SimulationResult)
    error("No implementation of union found for Simulation Results of type $(typeof(r1)) and $(typeof(r2))")
end
