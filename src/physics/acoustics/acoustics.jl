"""
    Acoustic{T<:AbstractFloat,Dim}(ρ::T, c::Complex{T})
    Acoustic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic acoustic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (1D) fields in Dim dimensions.
"""
struct Acoustic{T,Dim} <: PhysicalMedium{T,Dim,1}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

# basisorder_to_linearindices(::Type{Acoustic{T,3}}, order::Int) where T = (order^2 + 1):(order+1)^2
# basisorder_to_linearindices(::Type{Acoustic{T,2}}, order::Int) where T = 1:(2*order + 1)
basisorder_to_basislength(::Type{Acoustic{T,3}}, order::Int) where T = (order+1)^2
basisorder_to_basislength(::Type{Acoustic{T,2}}, order::Int) where T = 2*order + 1

basislength_to_basisorder(::Type{Acoustic{T,3}},len::Int) where T = Int(sqrt(len) - 1)
basislength_to_basisorder(::Type{Acoustic{T,2}},len::Int) where T = Int(T(len - 1) / T(2.0))

# Constructor which supplies the dimension without explicitly mentioning type
Acoustic(ρ::T,c::Union{T,Complex{T}},Dim::Integer) where {T} =  Acoustic{T,Dim}(ρ,Complex{T}(c))
Acoustic(Dim::Integer; ρ::T = 0.0, c::Union{T,Complex{T}} = 0.0) where {T} =  Acoustic{T,Dim}(ρ,Complex{T}(c))

import Base.show
function show(io::IO, p::Acoustic)
    # Acoustic template paramaters can be determined entirely from the medium and shape so we do not need to print them
    # Print is the style of the first constructor
    write(io, "Acoustic($(p.ρ), $(p.c), $(spatial_dimension(p)))")
    return
end

# Type aliases for convenience
TwoDimAcoustic{T} = Acoustic{T,2}

name(a::Acoustic{T,Dim}) where {Dim,T} = "$(Dim)D Acoustic"

"""
    impedance(medium::Acoustic)

Characteristic specific acoustic impedance (z₀) of medium
"""
impedance(medium::Acoustic) = medium.ρ * medium.c

function outgoing_basis_function(medium::Acoustic{T,2}, ω::T) where {T}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        k = ω/medium.c
        [hankelh1(m,k*r)*exp(im*θ*m) for m = -order:order]
    end
end

function outgoing_basis_function(medium::Acoustic{T,3}, ω::T) where {T}
    return function (order::Integer, x::AbstractVector{T})
        r, θ, φ  = cartesian_to_radial_coordinates(x)
        k = ω/medium.c

        Ys = spherical_harmonics(order, θ, φ)
        hs = [shankelh1(l,k*r) for l = 0:order]

        lm_to_n = lm_to_spherical_harmonic_index

        return [hs[l+1] * Ys[lm_to_n(l,m)] for l = 0:order for m = -l:l]
    end
end

function outgoing_translation_matrix(medium::Acoustic{T,2}, order::Integer, ω::T, x::AbstractVector{T}) where {T}
    translation_vec = outgoing_basis_function(medium, ω)(2order, x)
    N = basisorder_to_basislength(Acoustic{T,2},order)
    U = [translation_vec[n-m+N] for m in -order:order, n in -order:order]

    return U
end

function outgoing_translation_matrix(medium::Acoustic{T,3}, order::Integer, ω::T, x::AbstractVector{T}) where {T}
    us = outgoing_basis_function(medium, ω)(2*order,x)
    c = gaunt_coefficient

    ind(order::Int) = basisorder_to_basislength(Acoustic{T,3},order)
    U = [
        begin
            i1 = abs(l-dl) == 0 ? 1 : ind(abs(l-dl)-1) + 1
            i2 = ind(l+dl)

            cs = [c(T,l,m,dl,dm,l1,m1) for l1 = abs(l-dl):(l+dl) for m1 = -l1:l1]
            sum(us[i1:i2] .* cs)
        end
    for dl = 0:order for dm = -dl:dl for l = 0:order for m = -l:l];
    # U = [
    #     [(l,m),(dl,dm)]
    # for dl = 0:order for dm = -dl:dl for l = 0:order for m = -l:l]

    U = reshape(U, ((order+1)^2, (order+1)^2))

    return U
end

regular_basis_function(p::Particle{T,Dim,Acoustic{T,Dim}}, ω::T) where {T,Dim} = regular_basis_function(p.medium, ω)

function regular_basis_function(medium::Acoustic{T,2}, ω::Union{T,Complex{T}}) where T
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        k = ω/medium.c

        return [besselj(m,k*r)*exp(im*θ*m) for m = -order:order]
    end
end

function regular_basis_function(medium::Acoustic{T,3},  ω::Union{T,Complex{T}}) where T
    return function (order::Integer, x::AbstractVector{T})
        r, θ, φ  = cartesian_to_radial_coordinates(x)
        k = ω/medium.c

        Ys = spherical_harmonics(order, θ, φ)
        js = [sbesselj(l,k*r) for l = 0:order]

        lm_to_n = lm_to_spherical_harmonic_index

        return [js[l+1] * Ys[lm_to_n(l,m)] for l = 0:order for m = -l:l]
    end
end


# Check for material properties that don't make sense or haven't been implemented
"""
    check_material(p::Particle{T}, outer_medium::Acoustic{T})

Checks if wave scattering from the particle `p` is physically viable given the material properties of `p` and its surrounding medium `outer_medium`.
"""
function check_material(p::Particle{T}, outer_medium::Acoustic{T}) where T <: AbstractFloat

if isnan(abs(p.medium.c)*p.medium.ρ)
    throw(DomainError("Particle's phase speed times density is not a number!"))
elseif isnan(abs(outer_medium.c)*outer_medium.ρ)
    throw(DomainError("The medium's phase speed times density is not a number!"))
elseif iszero(outer_medium.c)
    throw(DomainError("Wave propagation in a medium with zero phase speed is not defined"))
elseif iszero(outer_medium.ρ) && iszero(p.medium.c*p.medium.ρ)
    throw(DomainError("Scattering in a medium with zero density from a particle with zero density or zero phase speed is not defined"))
elseif iszero(outer_radius(p))
    throw(DomainError("Scattering from a circle of zero radius is not implemented yet"))
end

return true

end

"""
    sound_hard([T::Type = Float64,] Dim::Integer)

Construct physical properties of a sound hard acoustic object with type T and dimension Dim.
Also known as [`rigid`](@ref) and equivalent to a [`zero_neumann`](@ref) pressure boundary condition.
"""
sound_hard(T::Type, Dim::Integer) = Acoustic{T,Dim}(T(Inf), one(T))

# If no type is given, assume Float64
sound_hard(Dim::Integer) = sound_hard(Float64, Dim)

"""
    hard(host_medium::Acoustic)

See [`sound_hard`](@ref).
"""
hard(host_medium::Acoustic{T,Dim}) where {T,Dim} = sound_hard(T, Dim)

"""
    rigid(host_medium::Acoustic)

See [`sound_hard`](@ref).
"""
rigid(host_medium::Acoustic{T,Dim}) where {T,Dim} = sound_hard(T, Dim)

"""
    zero_neumann(host_medium::Acoustic)

See [`sound_hard`](@ref).
"""
zero_neumann(host_medium::Acoustic{T,Dim}) where {T,Dim} = sound_hard(T, Dim)


"""
    sound_soft([T::Type = Float64,] Dim::Integer)

Construct physical properties of a sound hard acoustic object with type T and dimension Dim.
Equivalent to a [`zero_dirichlet`](@ref) pressure boundary condition.

"""
sound_soft(T::Type, Dim::Integer) = Acoustic{T,Dim}(zero(T), one(T))

# If no type is given, assume Float64
sound_soft(Dim::Integer) = sound_soft(Float64, Dim)

"""
    soft(host_medium::Acoustic)

See [`sound_soft`](@ref).
"""
soft(host_medium::Acoustic{T,Dim}) where {T,Dim} = sound_soft(T, Dim)

"""
    pressure_release(host_medium::Acoustic)

See [`sound_soft`](@ref).
"""
pressure_release(host_medium::Acoustic{T,Dim}) where {T,Dim} = sound_soft(T, Dim)

"""
    zero_dirichlet(host_medium::Acoustic)

See [`sound_soft`](@ref).
"""
zero_dirichlet(host_medium::Acoustic{T,Dim}) where {T,Dim} = sound_soft(T, Dim)
