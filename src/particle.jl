
struct Particle{T<:AbstractFloat,Dim,P<:PhysicalProperties,S<:Shape}
    medium::P
    shape::S
    # Enforce that the Dims and Types are all the same
    function Particle{T,Dim,P,S}(medium::P,shape::S) where {T,Dim,FieldDim,P<:PhysicalProperties{T,Dim,FieldDim},S<:Shape{T,Dim}}
        new{T,Dim,P,S}(medium,shape)
    end
end

# Shorthand for all Vectors of particles
Particles{T<:AbstractFloat,Dim} = Vector{Pt} where Pt<:(Particle{T,Dim,P,S} where S<:Shape where P<:PhysicalProperties)

# Convenience constructor which does not require explicit types/parameters
function Particle(medium::P,shape::S) where {Dim,FieldDim,T,P<:PhysicalProperties{T,Dim,FieldDim},S<:Shape{T,Dim}}
    Particle{T,Dim,P,S}(medium,shape)
end

# Shape hold infomation about origin of Particle
origin(p::Particle) = origin(p.shape)
boundary_points(p::Particle, num_points::Int = 3; kws...) = boundary_points(p.shape,num_points; kws...)

CircleParticle{T,P} = Particle{T,2,P,Circle{T}}

outer_radius(p::Particle) = outer_radius(p.shape)
volume(p::Particle) = volume(p.shape)

bounding_rectangle(p::Particle) = bounding_rectangle(p.shape)
bounding_rectangle(ps::Particles) = bounding_rectangle([p.shape for p in ps])

function volume(particles::Particles)
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
