"""
    Particle(medium::PhysicalMedium, shape::Shape)

Create particle with inner medium and shape (dimensions must agree).
"""
struct Particle{Dim,P<:PhysicalMedium,S<:Shape} <: AbstractParticle{Dim}
    medium::P
    shape::S
    # Enforce that the Dims are all the same
    function Particle{Dim,P,S}(medium::P,shape::S) where {Dim,FieldDim,P<:PhysicalMedium{Dim,FieldDim},S<:Shape{Dim}}
        new{Dim,P,S}(medium,shape)
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
struct CapsuleParticle{Dim,P<:PhysicalMedium,S<:Shape} <: AbstractParticle{Dim}
    outer::Particle{Dim,P,S}
    inner::Particle{Dim,P,S}
    # Enforce that particles are concentric
    function CapsuleParticle{Dim,P,S}(p2::Particle{Dim,P,S},p1::Particle{Dim,P,S}) where {Dim,P<:PhysicalMedium{Dim},S<:Shape}
        if origin(p1) != origin(p2) error("outer and inner particles should share the same origin") end
        if outer_radius(p1) >= outer_radius(p2)
            new{Dim,P,S}(p1,p2)
        else
            new{Dim,P,S}(p2,p1)
        end
    end
end

# Shorthand for all Vectors of particles
AbstractParticles{Dim} = Vector{Pt} where Pt<:AbstractParticle{Dim}

# Convenience constructor which does not require explicit types/parameters
function Particle(medium::P,s::S) where {Dim,P<:PhysicalMedium{Dim},S<:Shape{Dim}}
    Particle{Dim,P,S}(medium,s)
end

"""
    Particle(medium, radius)

Returns a particle shaped like a sphere or circle, when the particle shape is not given and with the specified `radius`.
"""
function Particle(medium::P, radius) where {Dim, P <: PhysicalMedium{Dim}}
    Particle{Dim,P,Sphere{typeof(radius),Dim}}(medium,Sphere(Dim,radius))
end

function CapsuleParticle(p1::Particle{Dim,P,S},p2::Particle{Dim,P,S}) where {Dim,S<:Shape,P<:PhysicalMedium}
    CapsuleParticle{Dim,P,S}(p1,p2)
end

shape(p::Particle) = p.shape
shape(p::CapsuleParticle) = p.outer.shape

# Shape holds infomation about origin of Particle
origin(p::AbstractParticle) = origin(shape(p))

boundary_points(p::AbstractParticle, num_points::Int = 3; kws...) = boundary_points(shape(p),num_points; kws...)

CircleParticle{T,P} = Particle{2,P,Sphere{T,2}}

outer_radius(p::AbstractParticle) = outer_radius(shape(p))
volume(p::AbstractParticle) = volume(shape(p))

bounding_box(p::AbstractParticle) = bounding_box(shape(p))
bounding_box(ps::AbstractParticles) = bounding_box([shape(p) for p in ps])

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

"""
    estimate_regular_basis_order(medium::P, ω::Number, radius::Number; tol = 1e-6)
"""
function estimate_outgoing_basisorder(medium::PhysicalMedium, p::Particle, ω::Number; tol = 1e-6)

    k = ω / real(medium.c)

    # A large initial guess
    L = Int(round(4 * abs(k*outer_radius(p)))) + 1

    ts = [
        norm(diag(t_matrix(p, medium, ω, l)))
    for l = 1:L]
    ts = ts ./ ts[1];

    l = findfirst(diff(ts) .< tol)

    return l
end

# Add generic documentation from shape
@doc (@doc issubset(::Shape, ::Shape)) issubset(::Shape, ::AbstractParticle)
@doc (@doc issubset(::Shape, ::Shape)) issubset(::AbstractParticle, ::Shape)
@doc (@doc issubset(::Shape, ::Shape)) issubset(::AbstractParticle, ::AbstractParticle)
