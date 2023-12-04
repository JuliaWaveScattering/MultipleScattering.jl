"""
Represent any source (incident) wave

Subtypes may have a symmetry (such as [`PlaneSource`](@ref)) and will contain information about physical medium.
"""
abstract type AbstractSource{P} end

"""
    PlaneSource(medium::P, amplitude::SVector, direction::SVector)

Is a struct type which describes a plane-wave source that drives/forces the whole system. It has three fields: a physical `medium`, an `amplitude` of the source, and the direction the propagate in `direction`.

For any given angular frequency ω, the PlaneSource has the value ``e^{i ω/c \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``c`` is the medium wavespeed and ``\\mathbf v`` is the direction.
"""
struct PlaneSource{T<:Real,Dim,FieldDim,P<:PhysicalMedium} <: AbstractSource{P}
    medium::P
    direction::SVector{Dim,T}
    position::SVector{Dim,T}
    amplitude::SVector{FieldDim,Complex{T}}
    # Check that P has same Dim and FieldDim
    function PlaneSource(medium::P, direction::AbstractArray, position = SVector(zeros(eltype(direction),Dim)...), amplitude = SVector{FieldDim}(ones(eltype(direction),FieldDim))) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim}}

        normw = sqrt(sum(direction .^2)) # note if direction is complex this is different from norm(direction)
        if !(normw ≈ 1)
            @warn "The direction will be normalised so that sum(direction .^2) == 1.0"
        end

        if length(direction) != Dim || length(position) != Dim
            throw(DimensionMismatch("The spatial dimensions of the medium do not match dimensions of the direction or position vector."))
        elseif length(amplitude) != FieldDim
            throw(DimensionMismatch("The dimensions of the field in the medium (i.e. always 1 for acoustics) do not match the dimension of the amplitude vector."))
        end

        T = promote_type(eltype(direction), eltype(position), real(eltype(amplitude)))
        new{T,Dim,FieldDim,P}(medium, SVector{Dim}(direction ./ normw), SVector{Dim}(position), complex(SVector{FieldDim}(amplitude)))
    end
end

field(s::PlaneSource{T}) where T = function (x::AbstractArray{T}, ω::T)
        s.amplitude * exp(im * (ω / s.medium.c) * dot(x - s.position,s.direction))
    end

field(s::PlaneSource{T,Dim,1}) where {T, Dim} = function (x::AbstractArray{T}, ω::T)
        s.amplitude[1] * exp(im * (ω / s.medium.c) * dot(x - s.position,s.direction))
    end

field(s::PlaneSource{T}, x::AbstractArray{T}, ω::T) where T = field(s)(x,ω)

function PlaneSource(medium::P;
        direction = vcat(SVector(1), SVector{Dim-1}(zeros(Float64,Dim-1))),
        position = zeros(eltype(direction),Dim),
        amplitude = SVector{FieldDim}(ones(eltype(direction),FieldDim))
    ) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim}}

    PlaneSource(medium,direction,position,amplitude)
end

Symmetry(s::PlaneSource{T,Dim}) where {Dim,T} = (abs(dot(s.direction,azimuthalnormal(Dim))) == norm(s.direction)) ? PlanarAzimuthalSymmetry{Dim}() : PlanarSymmetry{Dim}()

"""
    RegularSource(medium::P, field::Function, coef::Function)

Is a struct type which describes the source field that drives/forces the whole system. It is also known as an incident wave. It has three fields `RegularSource.medium`, `RegularSource.field`, and `RegularSource.coef`.

The source field at the position 'x' and angular frequency 'ω' is given by
```julia-repl
x = [1.0,0.0]
ω = 1.0
RegularSource.field(x,ω)
```

The field `RegularSource.coef`
regular_basis_function(medium::Acoustic{T,2}, ω::T)
"""
struct RegularSource{P<:PhysicalMedium,S<:AbstractSymmetry} <: AbstractSource{P}
    medium::P
    "Use: field(x,ω)"
    field::Function
    "Use: coefficients(n,x,ω)"
    coefficients::Function
    # Enforce that the Types are the same
    function RegularSource{P,S}(medium::P,field::Function,coefficients::Function) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim},S<:AbstractSymmetry{Dim}}
        s = new{P,S}(medium,field,coefficients)
        self_test(s)
        return s
    end
end

function RegularSource(medium::P,field::Function,coefficients::Function; symmetry::AbstractSymmetry{Dim} = WithoutSymmetry{Dim}()) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim}}
    return RegularSource{P,typeof(symmetry)}(medium,field,coefficients)
end

Symmetry(reg::RegularSource{P,S}) where {P,S} = S()

"""
    regular_spherical_coefficients(source::RegularSource)

return a function which can calculate the coefficients of a regular spherical wave basis.
"""
function regular_spherical_coefficients(source::RegularSource)
    source.coefficients
end

"""
Check that the source functions return the correct types
"""
function self_test(source::RegularSource{P}) where {P}
    ω = 1.0

    # choose rand postion, hopefully not the source position/origin
    x = SVector(rand(0.1:0.1:1.0,spatial_dimension(P))...)

    # Check that the result of field has the same dimension as field dimension of P
    if field_dimension(P) == 1
        length(first(source.coefficients(1,x,ω))) == field_dimension(P)
    # else # this else is not necessarily correct..
    #     source.coefficients(1,x,ω)::SVector{field_dimension(P),Complex{T}}
    end

    return true
end

field(s::RegularSource) = s.field
field(s::RegularSource, x::AbstractArray, ω::Number) = field(s)(x,ω)

function constant_source(medium::P, num::Number = zero(Float64) * im) where {P}
    return RegularSource(medium, (x,ω) -> num, (order,x,ω) -> [num])
end

"""
    regular_spherical_source(PhysicalMedium,regular_coefficients; amplitude = one(T), position = zeros(T,Dim))

``c_n = ```regular_coefficients[n]`, ``x_o=```position`, and let ``v_n(kx)`` represent the regular spherical basis with wavenumber ``k`` at position ``x``. The above returns a `source` where `source.field(x,ω) =```\\sum_n c_n v_n(kx -k x_0)``
"""
function regular_spherical_source(medium::PhysicalMedium{Dim},regular_coefficients::AbstractVector{CT};
        symmetry::AbstractSymmetry{Dim} = WithoutSymmetry{Dim}(),
        amplitude::Union{T,Complex{T}} = one(T),
        position::AbstractArray{T} = SVector(zeros(T,Dim)...)) where {T,Dim,CT<:Union{T,Complex{T}}}

    coeff_order = basislength_to_basisorder(typeof(medium),length(regular_coefficients))

    function source_field(x,ω)
        vs = regular_basis_function(medium, ω)(coeff_order,x-position)
        return amplitude * sum(vs .* regular_coefficients)

    end

    function source_coef(order,centre,ω)
        k = ω / medium.c

        # for the translation matrix below, order is the number columns and coeff_order is the number rows
        V = regular_translation_matrix(medium, order, coeff_order, ω, centre - position)

        return amplitude .* (transpose(V) * regular_coefficients)
    end

    return RegularSource(medium, source_field, source_coef; symmetry = symmetry)
end

"""
    source_expand(RegularSource, centre; basis_order = 4)

Returns a function of `(x,ω)` which approximates the value of the source at `(x,ω)`. That is, the source is written in terms of a regular basis expansion centred at `centre`.
"""
function source_expand(source::AbstractSource, centre::AbstractVector{T}; basis_order::Int = 4) where T

    # Convert to SVector for efficiency and consistency
    centre = SVector(centre...)

    return function (x::AbstractVector{T}, ω::T)
        vs = regular_basis_function(source.medium, ω)
        regular_coefficients = regular_spherical_coefficients(source)

        res = sum(regular_coefficients(basis_order,centre,ω) .* vs(basis_order, x - centre), dims = 1)[:]
        if length(res) == 1
            res = res[1]
        end
        
        return res
    end
end

import Base.(+)
function +(s1::RegularSource{P},s2::RegularSource{P})::RegularSource{P} where {P}
    if s1.medium != s2.medium
        error("Can not add sources from different physical mediums.")
    end
    field(x,ω) = s1.field(x,ω) + s2.field(x,ω)
    coef(n,centre,ω) = s1.coefficients(n,centre,ω) + s2.coefficients(n,centre,ω)

    sym1 = Symmetry(s1)
    sym2 = Symmetry(s2)
    S = typeof(Symmetry(sym1,sym2))

    RegularSource{P,S}(s1.medium,field,coef)
end

import Base.(*)
function *(a,s::RegularSource{P,S})::RegularSource{P,S} where {P,S}

    field(x,ω) = a * s.field(x,ω)
    coef(n,centre,ω) = a .* s.coefficients(n,centre,ω)

    RegularSource{P,S}(s.medium,field,coef)
end

*(s::RegularSource,a) = *(a,s)
