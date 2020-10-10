"""
    PhysicalMedium{T<:AbstractFloat,Dim,FieldDim}

An abstract type used to represent the physical medium, the dimension of
the field, and the number of spatial dimensions.
"""
abstract type PhysicalMedium{T<:AbstractFloat,Dim,FieldDim} end

"Extract the dimension of the field of this physical property"
field_dimension(p::PhysicalMedium{T,Dim,FieldDim}) where {T,Dim,FieldDim} = FieldDim

"Extract the dimension of the field of this type of physical property"
field_dimension(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalMedium{T,Dim,FieldDim}} = FieldDim

"Extract the dimension of the space that this physical property lives in"
spatial_dimension(p::PhysicalMedium{T,Dim,FieldDim}) where {Dim,FieldDim,T} = Dim

"Extract the dimension of the space that this type of physical property lives in"
spatial_dimension(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalMedium{T,Dim,FieldDim}} = Dim

"""
A basis for regular functions, that is, smooth functions. A series expansion in this basis should converge to any regular function within a ball.
"""
function regular_basis_function(medium::P, ω::T) where {P<:PhysicalMedium,T}
    error("No regular basis function implemented for this physics type.")
end

"""
Basis of outgoing wave. A series expansion in this basis should converge to any scattered field outside of a ball which contains the scatterer.
"""
function outgoing_basis_function(medium::P, ω::T) where {P<:PhysicalMedium,T}
    error("No outgoing basis function implmented for this physics type.")
end

"""
the field inside an AbstractParticle a some given point x.
"""
internal_field

"""
A tuples of vectors of the field close to the boundary of the shape. The field is calculated from sim::FrequencySimulation, but the PhysicalMedium inside and outside of the shape are assumed to be given by inside_medium and outside_medium.
"""
boundary_data


# estimate_regular_basisorder(medium::P, ka) where P<:PhysicalMedium = estimate_regular_basisorder(P, ka)

"""
    estimate_regular_basis_order(::Type{PhysicalMedium}, ka; tol = 1e-6)

where `ka = 2π * a / λ` is a ratio between a length `a` and a wavelength `λ`.
"""
function estimate_regular_basisorder(medium::P, ka; tol = 1e-6) where P<:PhysicalMedium

    vs = regular_basis_function(medium, medium.c)

    # A very large initial guess
    L = Int(round(10 * abs(ka)))

    kxs = ka .* rand(spatial_dimension(medium),10)

    l = nothing
    while isnothing(l)
        meanvs = mean(abs.(vs(L, kxs[:,i])) for i in axes(kxs,2))
        normvs = [norm(meanvs[basisorder_to_linearindices(P,i)]) for i = 1:L]
        l = findfirst(normvs .< tol)
        L = L + Int(round(abs(ka))) + 1
    end

    return l
end
