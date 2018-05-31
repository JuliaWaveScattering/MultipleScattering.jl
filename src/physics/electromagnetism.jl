
"""
Physical properties for a homogenous isotropic electromagnetic medium. Produces
a vector (??) field in arbitrary dimensions.
"""
type Electromagnetic{Dim,T} <: PhysicalProperties{T,Dim,3}
    μ::Complex{T} # Permeability
    ε::Complex{T} # Permittivity
    σ::Complex{T} # Conductivity
end

name(e::Electromagnetic{T,Dim}) where {T,Dim} = "$(Dim)D Electromagnetic"
