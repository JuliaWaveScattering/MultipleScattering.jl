
"Abstract class for Results of Simulations"
abstract type SimulationResult{T,Dim,FieldDim} end

"""
Struct to hold results of a FrequencySimulation
"""
struct FrequencySimulationResult{T<:AbstractFloat,Dim,FieldDim} <: SimulationResult{T,Dim,FieldDim}
    "Values of field through space (rows) and angular frequencies (columns)"
    field::Matrix{SVector{FieldDim,Complex{T}}}
    "Positions"
    x::Vector{SVector{Dim,T}}
    "Angular frequencies"
    ω::Vector{T}
end

function FrequencySimulationResult(field::Union{AbstractArray{Complex{T},2}, AbstractArray{AbstractVector{Complex{T}},2}}, x::AbstractVector{SVector{Dim,T}}, ω::AbstractVector{T}) where {Dim,T}

    field = [SVector(d...) for d in field]
    FieldDim = size(field[1],1)
    FrequencySimulationResult{T,Dim,FieldDim}(field, Vector(x), Vector(ω))
end

"""
Struct to hold results of a simulation through time
"""
struct TimeSimulationResult{T<:AbstractFloat,Dim,FieldDim} <: SimulationResult{T,Dim,FieldDim}
    "Values of field through space (rows) and time (columns)"
    field::Matrix{SVector{FieldDim,T}}
    "Positions"
    x::Vector{SVector{Dim,T}}
    "Times"
    t::Vector{T}
end

function TimeSimulationResult(time_field::Union{AbstractArray{T,2},AbstractArray{AbstractVector{T},2}}, x::AbstractVector{SVector{Dim,T}}, t::AbstractVector{T}) where {Dim,T}

    time_field = [SVector(d...) for d in time_field]
    FieldDim = size(time_field[1],1)
    TimeSimulationResult{T,Dim,FieldDim}(time_field, Vector(x), Vector(t))
end

"""
    field(result::SimulationResult, [i::Integer, j::Integer])

Get field from result, optionally specifying indices.

Returns single value of/matrix of complex SVectors() if vector field, and complex float if scalar field.
"""
field(result::SimulationResult) = result.field

function field(result::SimulationResult, i::Integer, j::Integer)
    result.field[i,j]
end

function field(result::SimulationResult{T,Dim,1}) where {Dim, T}
    map(x->x[1], result.field)
end

function field(result::SimulationResult{T,Dim,1}, i::Integer, j::Integer) where {Dim, T}
    result.field[i,j][1]
end

import Base.size
size(r::FrequencySimulationResult) = (length(r.x),length(r.ω))
size(r::TimeSimulationResult) = (length(r.x),length(r.t))

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

import Base.(+)
function +(s1::SimulationResult,s2::SimulationResult)::SimulationResult
    if typeof(s1) != typeof(s2)
        error("Can not sum different types $(typeof(s1)) and $(typeof(s2))")
    end

    return typeof(s1)(s1.field + s2.field, s1.x, getfield(s1,3))
end
function +(s::SimulationResult,a)::SimulationResult
    if typeof(s.field .+ a) != typeof(s.field)
        error("Summing SimulationResult with $a would cause SimulationResult.field to change its type.")
    end

    return typeof(s)(s.field .+ a, s.x, getfield(s,3))
end
+(a,s1::SimulationResult) = +(s1::SimulationResult,a)

import Base.(*)
function *(a,s::SimulationResult)::SimulationResult
    if typeof(s.field .* a) != typeof(s.field)
        error("Multiplying SimulationResult by $a would cause the field of SimulationResult to change type.")
    end

    return typeof(s)(s.field .* a, s.x, getfield(s,3))
end

*(s::SimulationResult,a) = *(a,s)
