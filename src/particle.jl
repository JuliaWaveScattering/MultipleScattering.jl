
type Particle{Dim,P<:PhysicalProperties,S<:Shape,T<:AbstractFloat}
    medium::P
    shape::S
    # Enforce that the Dims and Types are all the same
    function Particle{Dim,P,S,T}(medium::P,shape::S) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T},S<:Shape{Dim,T}}
        new{Dim,P,S,T}(medium,shape)
    end
end

# Convenience constructor which does not require explicit types/parameters
function Particle(medium::P,shape::S) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T},S<:Shape{Dim,T}}
    Particle{Dim,P,S,T}(medium,shape)
end

# Shape hold infomation about origin of Particle
origin(p::Particle) = origin(p.shape)
boundary_points(p::Particle, num_points::Int = 3; kws...) = boundary_points(p.shape,num_points; kws...)

CircleParticle{P, T} = Particle{2, P, Circle{T}, T}

outer_radius(p::Particle) = outer_radius(p.shape)
volume(p::Particle) = volume(p.shape)

function volume(particles::Vector{Pa}) where {Dim,Pa<:Particle{Dim}}
    mapreduce(volume, +, particles)
end

import Base.(==)
function ==(p1::Particle, p2::Particle)
    p1.medium == p2.medium &&
    p1.shape == p2.shape
end

"""
Retuns true if medium and shape of particles are the same, ignoring the origin
of shape
"""
function congruent(p1::Particle, p2::Particle)
    p1.medium == p2.medium &&
    congruent(p1.shape, p2.shape)
end


inside(shape::Shape, particle::Particle) = inside(shape, particle.shape)
