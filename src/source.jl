
type Source{P<:PhysicalProperties,T<:AbstractFloat}
    field::Function
    coef::Function
    # Enforce that the Types are the same
    function Source{P,T}(field::Function,coef::Function) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T}}
        s = new{P,T}(field,coef)
        self_test(s)
        return s
    end
end

"""
Returns a function that gives the value of the besselj expansion centred at centre
"""
function besselj_field(source::Source{Acoustic{2,T},T}, medium::Acoustic{2,T}, centre::AbstractVector{T}; hankel_order = 4) where T<:Number

    field(x,ω) = sum(
        source.coef(n,centre,ω)*besselj(n,ω/medium.c*norm(x - centre))*exp(im*n*atan2(x[2] - centre[2],x[1] - centre[1]))
    for n = -hankel_order:hankel_order)

   return field
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

function TwoDimAcousticPointSource{T}(medium::Acoustic{2,T}, source_position::AbstractVector{T}, amplitude::T = one(T))::Source{Acoustic{2,T},T}
    field(x,ω) = amplitude*Complex{T}(im/4)*hankelh1(0,ω/medium.c*norm(x-source_position))
    # using Graf's addition theorem
    coef(n,centre,ω) = amplitude*Complex{T}(im/4)*hankelh1(-n,ω/medium.c*norm(centre - source_position))*
exp(-Complex{T}(im)*n*atan2(centre[2] - source_position[2], centre[1] - source_position[1]))

    return Source{Acoustic{2,T},T}(field,coef)
end

function TwoDimAcousticPlanarSource{T}(medium::Acoustic{2,T}, source_position::AbstractVector{T}, source_direction::AbstractVector{T} = [one(T),zero(T)], amplitude::T = one(T))::Source{Acoustic{2,T},T}
    field(x,ω) = amplitude*exp(im*ω/medium.c*dot(x-source_position,source_direction))
    # Jacobi-Anger expansion
    coef(n,centre,ω) = field(centre,ω) * exp(im * n *(T(pi)/2 - atan2(source_direction[2],source_direction[1]) ))

    return Source{Acoustic{2,T},T}(field,coef)
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
