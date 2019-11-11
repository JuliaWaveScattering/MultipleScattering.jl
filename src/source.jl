
"""
    Source(medium::P, field::Function, coef::Function)

Is a struct type which describes the source field that drives/forces the whole system. It is also described as an incident wave. It has three fields `Source.medium`, `Source.field`, and `Source.coef`.

The source field at the position 'x' and angular frequency 'ω' is given by
```julia-repl
x = [1.0,0.0]
ω = 1.0
Source.field(x,ω)
```

The field `Source.coef`
regular_basis_function(medium::Acoustic{T,2}, ω::T)

"""
struct Source{P<:PhysicalProperties,T<:AbstractFloat}
    medium::P
    field::Function
    coef::Function
    # Enforce that the Types are the same
    function Source{P,T}(medium::P,field::Function,coef::Function) where {Dim,FieldDim,T,P<:PhysicalProperties{T,Dim,FieldDim}}
        s = new{P,T}(medium,field,coef)
        self_test(s)
        return s
    end
end

"""
Check that the source functions return the correct types
"""
function self_test(source::Source{P,T}) where {P,T}

    # Example data with correct dimensions and types from P and T
    x = SVector(ntuple(i->one(T),dim(P)))
    ω = one(T)

    # Check that the result of field has same dimension and type as PhysicalProperty field
    if field_dim(P) == 1
        source.field(x,ω)::Union{Complex{T},SVector{field_dim(P),Complex{T}}}
    else
        source.field(x,ω)::SVector{field_dim(P),Complex{T}}
    end

    # Check that the result of field has same dimension as field dimension of P
    if field_dim(P) == 1
        source.coef(1,x,ω)::Union{Vector{Complex{T}},Vector{SVector{field_dim(P),Complex{T}}}}
    # else # this else is not necessarily correct..
    #     source.coef(1,x,ω)::SVector{field_dim(P),Complex{T}}
    end

    return true
end

function constant_source(medium::P, num::Complex{T} = zero(Float64) * im) where {P,T}
    return Source{P,T}(medium, (x,ω) -> num, (order,x,ω) -> [num])
end

"""
    source_expand(Source, centre; basis_order = 4)

Returns a function of `(x,ω)` which approximates the value of the source at `(x,ω)`. That is, the source is written in terms of a regular basis expansion centred at `centre`.
"""
function source_expand(source::Source{P,T}, centre::AbstractVector{T}; basis_order::Int = 4) where {P,T}

    # Convert to SVector for efficiency and consistency
    centre = SVector(centre...)

    return function (x::AbstractVector{T}, ω::T)
        vs = regular_basis_function(source.medium, ω)
        sum(source.coef(basis_order,centre,ω) .* vs(basis_order, x - centre))
    end
end

import Base.(+)
function +(s1::Source{P,T},s2::Source{P,T})::Source{P,T} where {P,T}
    if s1.medium != s2.medium
        error("Can not add sources from different physical mediums.")
    end
    field(x,ω) = s1.field(x,ω) + s2.field(x,ω)
    coef(n,centre,ω) = s1.coef(n,centre,ω) + s2.coef(n,centre,ω)

    Source{P,T}(s1.medium,field,coef)
end

import Base.(*)
function *(a,s::Source{P,T})::Source{P,T} where {P,T}
    if typeof(one(Complex{T})*a) != Complex{T}
        error("Multiplying source by $a would cause source type to change, please explicitly cast $a to same type as Source ($T)")
    end

    field(x,ω) = a * s.field(x,ω)
    coef(n,centre,ω) = a .* s.coef(n,centre,ω)

    Source{P,T}(s.medium,field,coef)
end

*(s::Source,a) = *(a,s)
