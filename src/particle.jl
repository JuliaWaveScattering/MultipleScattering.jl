abstract type AbstractParticle{T,Dim} end

"A homogenous particle with any properties and shape"
struct Particle{T<:AbstractFloat,Dim,P<:PhysicalProperties,S<:Shape} <: AbstractParticle{T,Dim}
    medium::P
    shape::S
    # Enforce that the Dims and Types are all the same
    function Particle{T,Dim,P,S}(medium::P,shape::S) where {T,Dim,FieldDim,P<:PhysicalProperties{T,Dim,FieldDim},S<:Shape{T,Dim}}
        new{T,Dim,P,S}(medium,shape)
    end
end

"""
A particle within another particle, both with the same type of shape and origin. Produces a scalar (1D) field in arbitrary dimensions.
"""
struct CapsuleParticle{T<:AbstractFloat,Dim,P<:PhysicalProperties,S<:Shape} <: AbstractParticle{T,Dim}
    outer::Particle{T,Dim,P,S}
    inner::Particle{T,Dim,P,S}
    # Enforce that particles are concentric
    function CapsuleParticle{T,Dim,P,S}(p2::Particle{T,Dim,P,S},p1::Particle{T,Dim,P,S}) where {T,Dim,P<:PhysicalProperties{T,Dim},S<:Shape}
        if origin(p1) != origin(p2) error("outer and inner particles should share the same origin") end
        if outer_radius(p1) >= outer_radius(p2) error("inner particle needs to be smaller than outer particle") end
        new{T,Dim,P,S}(p1,p2)
    end
end

# Shorthand for all Vectors of particles
Particles{T<:AbstractFloat,Dim} = Vector{Pt} where Pt<:(Particle{T,Dim,P,S} where S<:Shape where P<:PhysicalProperties)

# Convenience constructor which does not require explicit types/parameters
function Particle(medium::P,shape::S) where {Dim,T,P<:PhysicalProperties{T,Dim},S<:Shape{T,Dim}}
    Particle{T,Dim,P,S}(medium,shape)
end

function CapsuleParticle(p1::Particle{T,Dim,P,S},p2::Particle{T,Dim,P,S}) where {T,Dim,S<:Shape,P<:PhysicalProperties}
    CapsuleParticle{T,Dim,P,S}(p1,p2)
end

shape(p::Particle) = p.shape
shape(p::CapsuleParticle) = p.outer.shape

# Shape holds infomation about origin of Particle
origin(p::AbstractParticle) = origin(shape(p))

boundary_points(p::AbstractParticle, num_points::Int = 3; kws...) = boundary_points(shape(p),num_points; kws...)

CircleParticle{T,P} = Particle{T,2,P,Circle{T}}

outer_radius(p::AbstractParticle) = outer_radius(shape(p))
volume(p::AbstractParticle) = volume(shape(p))

bounding_rectangle(p::AbstractParticle) = bounding_rectangle(shape(p))
bounding_rectangle(ps::Particles) = bounding_rectangle([shape(p) for p in ps])

function volume(particles::Particles)
    mapreduce(volume, +, particles)
end

import Base.(==)
function ==(p1::Particle, p2::Particle)
    p1.medium == p2.medium &&
    p1.shape == p2.shape
end

function ==(p1::CapsuleParticle, p2::CapsuleParticle)
    p1.outer == p2.outer &&
    p1.inner == p2.inner
end

import Base.isequal
function isequal(p1::Particle, p2::Particle)
    isequal(p1.medium, p2.medium) &&
    isequal(p1.shape, p2.shape)
end

"""
Retuns true if medium and shape of particles are the same, ignoring the origin
of shape
"""
function congruent(p1::Particle, p2::Particle)
    p1.medium == p2.medium &&
    congruent(p1.shape, p2.shape)
end

inside(shape::Shape, particle::AbstractParticle) = inside(shape, shape(particle))
