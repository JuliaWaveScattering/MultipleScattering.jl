"""
Object we can scatter waves off

Subtypes will contain information about shape and material properties. Most
crucially, they will implement the [`t_matrix`](@ref) function
"""
abstract type AbstractParticle{T,Dim} end

"""
    Particle(medium::PhysicalMedium, shape::Shape)

Create particle with inner medium and shape (types and dimension must agree).
"""
struct Particle{T<:AbstractFloat,Dim,P<:PhysicalMedium,S<:Shape} <: AbstractParticle{T,Dim}
    medium::P
    shape::S
    # Enforce that the Dims and Types are all the same
    function Particle{T,Dim,P,S}(medium::P,shape::S) where {T,Dim,FieldDim,P<:PhysicalMedium{T,Dim,FieldDim},S<:Shape{T,Dim}}
        new{T,Dim,P,S}(medium,shape)
    end
end

import Base.show
function show(io::IO, p::Particle)
    # Particle paramaters can be determined entirely from the medium and shape so we do not need to print them
    write(io, "Particle(")
    show(io, p.medium)
    write(io, ", ")
    show(io, p.shape)
    write(io, ")")
    return
end

"""
    CapsuleParticle(outer::Particle, inner::Particle)

A particle within another particle, both with the same shape type and origin.
"""
struct CapsuleParticle{T<:AbstractFloat,Dim,P<:PhysicalMedium,S<:Shape} <: AbstractParticle{T,Dim}
    outer::Particle{T,Dim,P,S}
    inner::Particle{T,Dim,P,S}
    # Enforce that particles are concentric
    function CapsuleParticle{T,Dim,P,S}(p2::Particle{T,Dim,P,S},p1::Particle{T,Dim,P,S}) where {T,Dim,P<:PhysicalMedium{T,Dim},S<:Shape}
        if origin(p1) != origin(p2) error("outer and inner particles should share the same origin") end
        if outer_radius(p1) >= outer_radius(p2)
            new{T,Dim,P,S}(p1,p2)
        else
            new{T,Dim,P,S}(p2,p1)
        end
    end
end

# Shorthand for all Vectors of particles
AbstractParticles{T<:AbstractFloat,Dim} = Vector{Pt} where Pt<:AbstractParticle{T,Dim}

# Convenience constructor which does not require explicit types/parameters
function Particle(medium::P,s::S) where {Dim,T,P<:PhysicalMedium{T,Dim},S<:Shape{T,Dim}}
    Particle{T,Dim,P,S}(medium,s)
end

function CapsuleParticle(p1::Particle{T,Dim,P,S},p2::Particle{T,Dim,P,S}) where {T,Dim,S<:Shape,P<:PhysicalMedium}
    CapsuleParticle{T,Dim,P,S}(p1,p2)
end

"Represents a set of particles."
struct Specie{T<:AbstractFloat,Dim,P<:AbstractParticle{T,Dim}}
    particle::P
    volume_fraction::T
    numberofparticles::Int
    exclusion_distance::T
end

# Convenience constructor which does not require explicit types/parameters
function Specie(p::AbstractParticle{T,Dim}; volume_fraction::T = 0.0, exclusion_distance::T = 1.005, numberofparticles::Int = -1) where {Dim,T<:AbstractFloat}
    Specie{T,Dim,typeof(p)}(p,volume_fraction,numberofparticles,exclusion_distance)
end

"Returns the volume fraction of the specie."
volume_fraction(s::Specie) = s.volume_fraction

"Returns the number density of the specie."
number_density(s::Specie{T,2,P}) where {T,P} = s.volume_fraction / (outer_radius(s.particle)^2 * pi)

number_density(s::Specie{T,3,P}) where {T,P} = s.volume_fraction / (T(4/3) * outer_radius(s.particle)^3 * pi)

shape(p::Particle) = p.shape
shape(p::CapsuleParticle) = p.outer.shape

# Shape holds infomation about origin of Particle
origin(p::AbstractParticle) = origin(shape(p))

boundary_points(p::AbstractParticle, num_points::Int = 3; kws...) = boundary_points(shape(p),num_points; kws...)

CircleParticle{T,P} = Particle{T,2,P,Circle{T}}

outer_radius(p::AbstractParticle) = outer_radius(shape(p))
volume(p::AbstractParticle) = volume(shape(p))

bounding_rectangle(p::AbstractParticle) = bounding_rectangle(shape(p))
bounding_rectangle(ps::AbstractParticles) = bounding_rectangle([shape(p) for p in ps])

function volume(particles::AbstractParticles)
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
    iscongruent(p1::AbstractParticle, p2::AbstractParticle)::Bool
    ≅(p1::AbstractParticle, p2::AbstractParticle)::Bool

Returns true if medium and shape of particles are the same, ignoring origin, false otherwise.
"""
iscongruent(p1::AbstractParticle, p2::AbstractParticle) = false # false by default, overload in specific examples

# Define synonym for iscongruent ≅, and add documentation
≅(p1::AbstractParticle, p2::AbstractParticle) = iscongruent(p1, p2)
@doc (@doc iscongruent(::AbstractParticle, ::AbstractParticle)) (≅(::AbstractParticle, ::AbstractParticle))

function iscongruent(p1::Particle, p2::Particle)
    p1.medium == p2.medium &&
    iscongruent(p1.shape, p2.shape)
end

function iscongruent(p1::CapsuleParticle, p2::CapsuleParticle)
    iscongruent(p1.inner, p2.inner) && iscongruent(p1.outer, p2.outer)
end

import Base.in
"""
    in(vector, particle)::Bool

Returns true if vector is in interior of shape of particle, false otherwise.
"""
in(x::AbstractVector, particle::AbstractParticle) = in(x, shape(particle))

import Base.issubset
issubset(s::Shape, particle::AbstractParticle) = issubset(s, shape(particle))
issubset(particle::AbstractParticle, s::Shape) = issubset(shape(particle), s)
issubset(particle1::AbstractParticle, particle2::AbstractParticle) = issubset(shape(particle1), shape(particle2))

# Add generic documentation from shape
@doc (@doc issubset(::Shape, ::Shape)) issubset(::Shape, ::AbstractParticle)
@doc (@doc issubset(::Shape, ::Shape)) issubset(::AbstractParticle, ::Shape)
@doc (@doc issubset(::Shape, ::Shape)) issubset(::AbstractParticle, ::AbstractParticle)
