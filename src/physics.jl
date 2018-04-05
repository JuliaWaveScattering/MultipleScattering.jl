"""
Holds information about the physical properties of the medium, the dimension of
the field and the number of dimensions it is a function of.
"""
abstract type PhysicalProperties{Dim,FieldDim,T<:AbstractFloat} end

"""
Physical properties for a homogenous isotropic acoustic medium. Produces a
scalar (1D) field in arbitrary dimensions.
"""
type Acoustics{Dim,T} <: PhysicalProperties{Dim,1,T}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

# Constructor which supplies the dimension without explicitly mentioning type
Acoustics(ρ::T,c::Complex{T},Dim::Integer) where {T} =  Acoustics{Dim,T}(c,ρ)

name(a::Acoustics{Dim,T}) where {Dim,T} = "$(Dim)D Acoustic"

function basis_function(medium::Acoustics{2,T}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::MVector{2,T})
        r = norm(x)
        θ = atan2(x[2],x[1])
        k = ω/medium.c
        hankelh1(k*r,m)*exp(im*θ*m)
    end
end

"""
Basis functions in a specific dimension for a specific physics type.
"""
function basis_function(medium::PhysicalProperties, ω::T) where {T}
    error("Basis functions not implmented for this physics type.")
end

"""
Physical properties for a capsule with the same inner and outer shape, where the
characteristic length scale of the two shapes is `size_ratio`. Inner and outer
are both a homogenous isotropic acoustic medium. Produces a scalar (1D) field in
arbitrary dimensions.
"""
type AcousticCapsule{T,Dim} <: PhysicalProperties{Dim,1,T}
    inner_ρ::T # Density in inner shape
    inner_c::Complex{T} # Phase velocity in inner shape
    outer_ρ::T # Density in outer shape
    outer_c::Complex{T} # Phase velocity in outer shape
    size_ratio::T # Ratio of size of inner shape to outer shape
end

name(a::AcousticCapsule{T,Dim}) where {Dim,T} = "$(Dim)D Acoustic Capsule"


"""
Physical properties for a homogenous isotropic electromagnetic medium. Produces
a vector (??) field in arbitrary dimensions.
"""
type Electromagnetism{Dim,T} <: PhysicalProperties{Dim,3,T}
    μ::Complex{T} # Permeability
    ε::Complex{T} # Permittivity
    σ::Complex{T} # Conductivity
end

name(e::Electromagnetism{T,Dim}) where {T,Dim} = "$(Dim)D Electromagnetic"
