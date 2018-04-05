
"Abstract class for Results of Simulations"
abstract type Result end

struct FrequencySimulationResult{P,T,Dim} <: Result where P <: PhysicalProperties{Dim,FieldDim,T}, T <: AbstractFloat, Dim::Int
    field::Matrix{MVector{T,FieldDim}}
    x::Vector{MVector{T,FieldDim}}
    ω::RowVector{T}
end

struct TimeSimulationResult{P,T,Dim} <: Result where P <: PhysicalProperties{Dim,FieldDim,T}, T <: AbstractFloat, Dim::Int
    field::Vector{Vector{{MVector{T,FieldDim}}}}
    x::Vector{MVector{Dim,T}}
    t::RowVector{T}
end

import Base.union

"""
Combine two FrequencyResults intelligently, allows user to optionally sort ω
and x
"""
function union(r1::FrequencySimulationResult, r2::FrequencySimulationResult; sort::function=identity)
    if r1.x == r2.x
    elseif r1.k == r2.k
    end
end

function union(r1::R1,r2::R2) where R1, R2 <: Result
    error("No implementation of union found for Results of type $(typeof(r1)) and $(typeof(r2))")
end
