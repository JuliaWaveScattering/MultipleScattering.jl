
abstract PhysicalProperties{T,Dim,FieldDim} end

type Acoustics{T,Dim} <: PhysicalProperties{T,Dim,1}
    ρ::T #the medium's density
    c::Complex{T} #the medium's phase velocity
end

function transmission_coefs(medium::Acoustics{T,Dim}, particle::Acoustics{T,Dim}, ω::T) where T <: AbstractFloat, Dim::Int
    return #Zns???
end

# Could this type of framework with multiple dispacth be used in full FSI??
function transmission_coefs(medium::Acoustics, particle::TimeHarmonicLinearElast)

end

function φ(x::MVector{2,T}, physics::Acoustics{T,2}, ω::T, m::Int)
    r = x[1]; θ = x[2]
    hankelh1(ω/physics.c*r,m)*exp(im*θ*m)
end

type Electromagnetism{T,Dim} <: PhysicalProperties{T,Dim,1}
    μ::Complex{T} # Permeability
    ε::Complex{T} # Permittivity
    σ::Complex{T} #Conductivity
end
