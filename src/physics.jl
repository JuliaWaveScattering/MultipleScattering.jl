"""
Holds information about the physical properties of the medium, the dimension of
the field and the number of dimensions it is a function of.
"""
abstract PhysicalProperties{Dim,FieldDim,T} where T <: AbstractFloat, Dim::Integer, FieldDim::Integer end

"""
Physical properties for a homogenous isotropic acoustic medium. Produces a
scalar (1D) field in arbitrary dimensions.
"""
type Acoustics{T,Dim} <: PhysicalProperties{T,Dim,1}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

function φ(x::MVector{2,T}, physics::Acoustics{T,2}, ω::T, m::Integer)
    r = norm(x)
    θ = atan2(x[2],x[1])
    k = ω/physics.c
    hankelh1(k*r,m)*exp(im*θ*m)
end

"""
Basis functions in a specific dimension for a specific physics type.
"""
function φ(x::MVector{Dim,T}, physics::P, ω::T, m::Int) where Dim::Integer, P <: PhysicalProperties{Dim,FieldDim,T}
    error("Basis functions not implmented for this physics type.")
end

"""
Physical properties for a capsule with the same inner and outer shape, where the
characteristic length scale of the two shapes is `size_ratio`. Inner and outer
are both a homogenous isotropic acoustic medium. Produces a scalar (1D) field in
arbitrary dimensions.
"""
type CapsuleAcoustics{T,Dim} <: PhysicalProperties{T,Dim,1}
    inner_ρ::T # Density in inner shape
    inner_c::Complex{T} # Phase velocity in inner shape
    outer_ρ::T # Density in outer shape
    outer_c::Complex{T} # Phase velocity in outer shape
    size_ratio::T # Ratio of size of inner shape to outer shape
end

"""
Physical properties for a homogenous isotropic electromagnetic medium. Produces
a vector (??) field in arbitrary dimensions.
"""
type Electromagnetism{T,Dim} <: PhysicalProperties{T,Dim,3}
    μ::Complex{T} # Permeability
    ε::Complex{T} # Permittivity
    σ::Complex{T} # Conductivity
end
