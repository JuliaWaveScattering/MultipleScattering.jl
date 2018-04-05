
type Particle{Dim,P<:PhysicalProperties,S<:Shape,T<:AbstractFloat}
    position::MVector{Dim,T}
    medium::P
    shape::S
    # Enforce that the Dims and Types are all the same
    function Particle{Dim,P,S,T}(position::MVector{Dim,T},medium::P,shape::S) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T},S<:Shape{Dim,T}}
        new{Dim,P,S,T}(position,medium,shape)
    end
end

# Convenience constructor which does not require explicit types/parameters
function Particle(position::MVector{Dim,T},medium::P,shape::S) where {Dim,T,P,S}
    Particle{Dim,P,S,T}(position,medium,shape)
end

# CircleParticle{P, T, Dim} = Particle{T, Dim, P, Circle{T}}

volume(p::Particle) = volume(p.shape)

volume(particles::Vector{Particle}) = mapreduce(volume, +, particles)

import Base.(==)
function ==(p1::Particle, p2::Particle)
    p1.position == p2.position &&
    p1.medium == p2.medium &&
    p1.shape == p2.shape
end
