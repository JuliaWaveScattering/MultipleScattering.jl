"""
    Acoustic{T<:AbstractFloat,Dim}(ρ::T, c::Complex{T})
    Acoustic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic acoustic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (1D) fields in Dim dimensions.
"""
struct Acoustic{T,Dim} <: PhysicalMedium{Dim,1}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

# Constructor which supplies the dimension without explicitly mentioning type
Acoustic(ρ::T,c::Union{T,Complex{T}},Dim::Integer) where {T<:Number} =  Acoustic{T,Dim}(ρ,Complex{T}(c))
Acoustic(Dim::Integer; ρ::T = 0.0, c::Union{T,Complex{T}} = 0.0) where {T<:Number} =  Acoustic{T,Dim}(ρ,Complex{T}(c))

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

# Check for material properties that don't make sense or haven't been implemented
"""
    check_material(p::Particle, outer_medium::Acoustic)

Checks if wave scattering from the particle `p` is physically viable given the material properties of `p` and its surrounding medium `outer_medium`.
"""
function check_material(p::Particle, outer_medium::Acoustic)

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
hard(host_medium::Acoustic{T,Dim}) where {Dim,T} = sound_hard(T, Dim)

"""
    rigid(host_medium::Acoustic)

See [`sound_hard`](@ref).
"""
rigid(host_medium::Acoustic{T,Dim}) where {Dim,T} = sound_hard(T, Dim)

"""
    zero_neumann(host_medium::Acoustic)

See [`sound_hard`](@ref).
"""
zero_neumann(host_medium::Acoustic{T,Dim}) where {Dim,T} = sound_hard(T, Dim)


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
soft(host_medium::Acoustic{T,Dim}) where {Dim,T} = sound_soft(T, Dim)

"""
    pressure_release(host_medium::Acoustic)

See [`sound_soft`](@ref).
"""
pressure_release(host_medium::Acoustic{T,Dim}) where {Dim,T} = sound_soft(T, Dim)

"""
    zero_dirichlet(host_medium::Acoustic)

See [`sound_soft`](@ref).
"""
zero_dirichlet(host_medium::Acoustic{T,Dim}) where {Dim,T} = sound_soft(T, Dim)

"""
    internal_field(x::AbstractVector, p::Particle{Dim,Acoustic{T,Dim}},  source::RegularSource, ω::T, scattering_coefficients::AbstractVector{Complex{T}})

The internal field for an acoustic particle in an acoustic medium. For a sphere and circular cylinder the result is exact, for everything else it is an approximation which assumes smooth fields.
"""
function internal_field(x::AbstractVector{T}, p::Particle{Dim,Acoustic{T,Dim}}, source::RegularSource{Acoustic{T,Dim}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}}) where {Dim,T}
    if !(x ∈ p)
        @error "Point $x is not inside the particle with shape $(p.shape)"
    end
    if iszero(p.medium.c) || isinf(abs(p.medium.c))
        return zero(Complex{T})
    else
        fs = scattering_coefficients
        order = basislength_to_basisorder(Acoustic{T,Dim},length(fs))
        r = outer_radius(p)

        t_mat = t_matrix(p, source.medium, ω, order)
        vs = regular_radial_basis(source.medium, ω, order, r)
        vos = regular_radial_basis(p.medium, ω, order, r)
        us = outgoing_radial_basis(source.medium, ω, order, r)

        internal_coefs = (vs .* (inv(t_mat) * fs) + us .* fs) ./ vos
        inner_basis = regular_basis_function(p, ω)

        return sum(inner_basis(order, x-origin(p)) .* internal_coefs)
    end
end
