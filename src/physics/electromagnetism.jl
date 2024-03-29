
"""
    Electromagnetic

Represents the physical properties for a homogenous isotropic electromagnetic medium. Produces
a three dimensional vector field.
"""
struct Electromagnetic{Dim,T} <: PhysicalMedium{Dim,3}
    μ::Complex{T} # Permeability
    ε::Complex{T} # Permittivity
    σ::Complex{T} # Conductivity
end

name(e::Electromagnetic{Dim,T}) where {Dim,T} = "$(Dim)D Electromagnetic"
