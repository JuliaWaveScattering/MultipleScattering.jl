
"""
Physical properties for a homogenous isotropic electromagnetic medium. Produces
a three dimensional vector field.
"""
struct Electromagnetic{Dim,T} <: PhysicalProperties{T,Dim,3}
    μ::Complex{T} # Permeability
    ε::Complex{T} # Permittivity
    σ::Complex{T} # Conductivity
end

name(e::Electromagnetic{T,Dim}) where {T,Dim} = "$(Dim)D Electromagnetic"
