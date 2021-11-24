"""
    PhysicalMedium{Dim,FieldDim}

An abstract type used to represent the physical medium, the dimension of
the field, and the number of spatial dimensions.
"""
abstract type PhysicalMedium{Dim,FieldDim} end

"Extract the dimension of the field of this physical property"
field_dimension(p::PhysicalMedium{Dim,FieldDim}) where {Dim,FieldDim} = FieldDim

"Extract the dimension of the field of this type of physical property"
field_dimension(p::Type{P}) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim}} = FieldDim

"Extract the dimension of the space that this physical property lives in"
spatial_dimension(p::PhysicalMedium{Dim,FieldDim}) where {Dim,FieldDim} = Dim

"Extract the dimension of the space that this type of physical property lives in"
spatial_dimension(p::Type{P}) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim}} = Dim

"""
A basis for regular functions, that is, smooth functions. A series expansion in this basis should converge to any regular function within a ball.
"""
regular_basis_function

"""
Basis of outgoing wave. A series expansion in this basis should converge to any scattered field outside of a ball which contains the scatterer.
"""
outgoing_basis_function

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
    estimate_regular_basis_order(medium::PhysicalMedium, ω::Number, radius::Number; tol = 1e-6)
"""
function estimate_regular_basisorder(medium::PhysicalMedium, ω::Number, radius::Number; tol = 1e-6)

    @error "This is not complete"

    k = ω / real(medium.c)
    vs = regular_basis_function(medium, ω)

    # A large initial guess
    L = Int(round(4 * abs(k*radius)))

    xs = radius .* rand(spatial_dimension(medium),10)

    l = nothing
    while isnothing(l)
        meanvs = [
            mean(norm(vs(l, xs[:,i])) for i in axes(xs,2))
        for l = 1:L]

        meanvs = mean(abs.(vs(L, xs[:,i])) for i in axes(xs,2))
        normvs = [
            norm(meanvs[basisorder_to_basislength(P,i-1):basisorder_to_basislength(P,i)])
        for i = 1:L]
        l = findfirst(normvs .< tol)
        L = L + Int(round(abs(k * radius / 2.0))) + 1
    end

    return l
end
