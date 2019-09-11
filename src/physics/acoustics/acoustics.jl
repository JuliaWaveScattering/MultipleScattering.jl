"""
    Acoustic{T<:AbstractFloat,Dim}(ρ::T, c::Complex{T})
    Acoustic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic acoustic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (1D) fields in Dim dimensions.
"""
struct Acoustic{T,Dim} <: PhysicalProperties{T,Dim,1}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

# Constructor which supplies the dimension without explicitly mentioning type
Acoustic(ρ::T,c::Complex{T},Dim::Integer) where {T} =  Acoustic{T,Dim}(ρ,c)
Acoustic(ρ::T,c::T,Dim::Integer) where {T} =  Acoustic{T,Dim}(ρ,Complex{T}(c))
Acoustic(Dim::Integer; ρ::T = 0.0, c::T = 0.0) where {T} =  Acoustic{T,Dim}(ρ,Complex{T}(c))

import Base.show
function show(io::IO, p::Acoustic)
    # Acoustic template paramaters can be determined entirely from the medium and shape so we do not need to print them
    # Print is the style of the first constructor
    write(io, "Acoustic($(p.ρ), $(p.c), $(dim(p)))")
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

include("circle.jl")
include("concentric_capsule.jl")
include("source.jl")
include("boundary_data.jl")


function outgoing_basis_function(medium::Acoustic{T,2}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::SVector{2,T})
        r = norm(x)
        θ = atan(x[2],x[1])
        k = ω/medium.c
        hankelh1(m,k*r)*exp(im*θ*m)
    end
end

"Basis function when inside a particle. Assumes particle is a circle, which approximately works for all shapes."
function regular_basis_function(p::Particle{T,2,Acoustic{T,2}}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::SVector{2,T})
        r = norm(x)
        θ = atan(x[2],x[1])
        k = ω/p.medium.c
        besselj(m,k*r)*exp(im*θ*m)
    end
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
