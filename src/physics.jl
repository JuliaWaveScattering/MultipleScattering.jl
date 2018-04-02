"""
Holds information about the physical properties of the medium, the dimension of
the field and the number of dimensions it is a function of.
"""
abstract PhysicalProperties{T,Dim,FieldDim} where T <: AbstractFloat, Dim::Integer, FieldDim::Integer end

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
function φ(x::MVector{Dim,T}, physics::P, ω::T, m::Int) where Dim::Integer, P <: PhysicalProperties{T,Dim,FieldDim}
    error("Basis functions not implmented for this physics type.")
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
