"""
    PhysicalMedium{Dim,FieldDim}

An abstract type used to represent the physical medium, the dimension of
the field, and the number of spatial dimensions. We expect every sub-type of PhysicalMedium{Dim,1} to have a field c which is the complex wave speed. 
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
    ScalarMedium{Dim}

An type used to represent a scalar wave which satisfies a Helmoholtz equations. Dim represents the number of spatial dimensions.
"""
struct ScalarMedium{T,Dim} <: PhysicalMedium{Dim,1}
    c::Complex{T} # complex wavespeed
end

"""
A basis for regular functions, that is, smooth functions. A series expansion in this basis should converge to any regular function within a ball.
"""
regular_basis_function

"""
Basis of outgoing wave when the dependance on the angles has been removed.
"""
outgoing_radial_function

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

## Below are many of the essential functions for scalar wave potentials (i.e. PhysicalMedium{Dim,1}) which are used to build all other types of vector waves (as we strict ourselves to isotropy.)

# basisorder_to_linearindices(::Type{Acoustic{T,3}}, order::Int) where T = (order^2 + 1):(order+1)^2
# basisorder_to_linearindices(::Type{Acoustic{T,2}}, order::Int) where T = 1:(2*order + 1)
basisorder_to_basislength(::Type{P}, order::Int) where P <: PhysicalMedium{3,1} = (order+1)^2
basisorder_to_basislength(::Type{P}, order::Int) where P <: PhysicalMedium{2,1} = 2*order + 1

basislength_to_basisorder(::Type{P},len::Int) where P <: PhysicalMedium{3,1} = Int(sqrt(len) - 1)
basislength_to_basisorder(::Type{P},len::Int) where P <: PhysicalMedium{2,1} = Int((len - 1) / 2.0)

function outgoing_radial_basis(medium::PhysicalMedium{2,1}, ω::T, order::Integer, r::T) where {T<:Number}
    k = ω/medium.c
    return transpose(hankelh1.(-order:order,k*r))
end

function outgoing_basis_function(medium::PhysicalMedium{2,1}, ω::T) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        k = ω/medium.c
        [hankelh1(m,k*r)*exp(im*θ*m) for m = -order:order] |> transpose
    end
end

function outgoing_radial_basis(medium::PhysicalMedium{3,1}, ω::T, order::Integer, r::T) where {T<:Number}
    k = ω/medium.c
    hs = shankelh1.(0:order,k*r)
    return  [hs[l+1] for l = 0:order for m = -l:l] |> transpose
end

function outgoing_basis_function(medium::PhysicalMedium{3,1}, ω::T) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ, φ  = cartesian_to_radial_coordinates(x)
        k = ω/medium.c

        Ys = spherical_harmonics(order, θ, φ)
        hs = [shankelh1(l,k*r) for l = 0:order]

        lm_to_n = lm_to_spherical_harmonic_index

        return [hs[l+1] * Ys[lm_to_n(l,m)] for l = 0:order for m = -l:l] |> transpose
    end
end

function outgoing_translation_matrix(medium::PhysicalMedium{2,1}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}
    translation_vec = outgoing_basis_function(medium, ω)(in_order + out_order, x)
    U = [
        translation_vec[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    return U
end

function regular_translation_matrix(medium::PhysicalMedium{2,1}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}
    translation_vec = regular_basis_function(medium, ω)(in_order + out_order, x)
    V = [
        translation_vec[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    return V
end

function outgoing_translation_matrix(medium::PhysicalMedium{3,1}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}

    us = outgoing_basis_function(medium, ω)(in_order + out_order,x)
    c = gaunt_coefficient

    ind(order::Int) = basisorder_to_basislength(Acoustic{T,3},order)
    U = [
        begin
            i1 = abs(l-dl) == 0 ? 1 : ind(abs(l-dl)-1) + 1
            i2 = ind(l+dl)

            cs = [c(T,l,m,dl,dm,l1,m1) for l1 = abs(l-dl):(l+dl) for m1 = -l1:l1]
            sum(us[i1:i2] .* cs)
        end
    for dl = 0:in_order for dm = -dl:dl for l = 0:out_order for m = -l:l];
    # U = [
    #     [(l,m),(dl,dm)]
    # for dl = 0:order for dm = -dl:dl for l = 0:order for m = -l:l]

    U = reshape(U, ((out_order+1)^2, (in_order+1)^2))

    return U
end

# NOTE that medium in both functions below is only used to get c and to identify typeof(medium)
regular_basis_function(medium::PhysicalMedium{Dim,1}, ω::Union{T,Complex{T}}) where {Dim,T} = regular_basis_function(ω/medium.c, medium)

function regular_radial_basis(medium::PhysicalMedium{2,1}, ω::T, order::Integer, r::T) where {T<:Number}
    k = ω / medium.c
    return transpose(besselj.(-order:order,k*r))
end

function regular_basis_function(wavenumber::Union{T,Complex{T}}, ::PhysicalMedium{2,1}) where T
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        k = wavenumber

        return [besselj(m,k*r)*exp(im*θ*m) for m = -order:order] |> transpose
    end
end

function regular_radial_basis(medium::PhysicalMedium{3,1}, ω::T, order::Integer, r::T) where {T<:Number}
    k = ω / medium.c
    js = sbesselj.(0:order,k*r)

    return [js[l+1] for l = 0:order for m = -l:l] |> transpose
end

function regular_basis_function(wavenumber::Union{T,Complex{T}}, ::PhysicalMedium{3,1}) where T
    return function (order::Integer, x::AbstractVector{T})
        r, θ, φ  = cartesian_to_radial_coordinates(x)

        Ys = spherical_harmonics(order, θ, φ)
        js = [sbesselj(l,wavenumber*r) for l = 0:order]

        lm_to_n = lm_to_spherical_harmonic_index

        return [js[l+1] * Ys[lm_to_n(l,m)] for l = 0:order for m = -l:l] |> transpose
    end
end

function regular_translation_matrix(medium::PhysicalMedium{3,1}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}
    vs = regular_basis_function(medium, ω)(in_order + out_order,x)
    c = gaunt_coefficient

    ind(order::Int) = basisorder_to_basislength(Acoustic{T,3},order)
    V = [
        begin
            i1 = (abs(l-dl) == 0) ? 1 : ind(abs(l-dl)-1) + 1
            i2 = ind(l+dl)

            cs = [c(T,l,m,dl,dm,l1,m1) for l1 = abs(l-dl):(l+dl) for m1 = -l1:l1]
            sum(vs[i1:i2] .* cs)
        end
    for dl = 0:in_order for dm = -dl:dl for l = 0:out_order for m = -l:l];

    V = reshape(V, ((out_order+1)^2, (in_order+1)^2))

    return V
end


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
