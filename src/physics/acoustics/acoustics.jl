"""
Physical properties for a homogenous isotropic acoustic medium. Produces a
scalar (1D) field in arbitrary dimensions.
"""
struct Acoustic{T,Dim} <: PhysicalProperties{T,Dim,1}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

# Constructor which supplies the dimension without explicitly mentioning type
Acoustic(ρ::T,c::Complex{T},Dim::Integer) where {T} =  Acoustic{T,Dim}(ρ,c)
Acoustic(ρ::T,c::T,Dim::Integer) where {T} =  Acoustic{T,Dim}(ρ,Complex{T}(c))

# Type aliases for convenience
TwoDimAcoustic{T} = Acoustic{T,2}

name(a::Acoustic{T,Dim}) where {Dim,T} = "$(Dim)D Acoustic"

include("circle.jl")
include("capsule.jl")
include("source.jl")
include("boundary_data.jl")


function basis_function(medium::Acoustic{T,2}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::SVector{2,T})
        r = norm(x)
        θ = atan2(x[2],x[1])
        k = ω/medium.c
        hankelh1(m,k*r)*exp(im*θ*m)
    end
end

"Basis function when inside a particle. Assumes particle is a circle, which approximately works for all shapes."
function basis_function(p::Particle{T,2,Acoustic{T,2}}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::SVector{2,T})
        r = norm(x)
        θ = atan2(x[2],x[1])
        k = ω/p.medium.c
        besselj(m,k*r)*exp(im*θ*m)
    end
end
