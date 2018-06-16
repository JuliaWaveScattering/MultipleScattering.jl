
"""
    Source{P,T}(field::Function,coef::Function) where P<:PhysicalProperties{T,Dim,FieldDim}

A struct which describes the source that drives/forces the whole system. It is also described as an incident wave.

```julia-repl
x = [1.0,0.0]
ω = 1.0
Source.field(x,ω)
```
should return the wave field at position 'x' and angular frequency 'ω'
"""
struct Source{P<:PhysicalProperties,T<:AbstractFloat}
    field::Function
    coef::Function
    # Enforce that the Types are the same
    function Source{P,T}(field::Function,coef::Function) where {Dim,FieldDim,T,P<:PhysicalProperties{T,Dim,FieldDim}}
        s = new{P,T}(field,coef)
        self_test(s)
        return s
    end
end

"""
Check that the source functions return the correct types
"""
function self_test(source::Source{P,T}) where {P,T}

    # Example data with correct dimensions and types from P and T
    x = SVector(ntuple(i->zero(T),dim(P)))
    ω = one(T)

    # Check that the result of field has same dimension and type as PhysicalProperty field
    if field_dim(P) == 1
        source.field(x,ω)::Union{Complex{T},SVector{field_dim(P),Complex{T}}}
    else
        source.field(x,ω)::SVector{field_dim(P),Complex{T}}
    end

    # Check that the result of field has same dimension as field dimension of P
    if field_dim(P) == 1
        source.coef(1,x,ω)::Union{Complex{T},SVector{field_dim(P),Complex{T}}}
    else
        source.coef(1,x,ω)::SVector{field_dim(P),Complex{T}}
    end

    return true
end

import Base.(+)
function +(s1::Source{P,T},s2::Source{P,T})::Source{P,T} where {P,T}
    field(x,ω) = s1.field(x,ω) + s2.field(x,ω)
    coef(n,centre,ω) = s1.coef(n,centre,ω) + s2.coef(n,centre,ω)

    Source{P,T}(field,coef)
end

import Base.(*)
function *(a,s::Source{P,T})::Source{P,T} where {P,T}
    if typeof(one(Complex{T})*a) != Complex{T}
        error("Multiplying source by $a would cause source type to change, please explicitly cast $a to same type as Source ($T)")
    end

    field(x,ω) = a*s.field(x,ω)
    coef(n,centre,ω) = a*s.coef(n,centre,ω)

    Source{P,T}(field,coef)
end

*(s::Source,a) = *(a,s)
