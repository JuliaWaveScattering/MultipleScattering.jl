
type Particle{T,Dim,P,S} where P <: PhysicalProperties{Dim,FieldDim,T}, S <: Shape{Dim, T}
    position::MVector{Dim,T}
    medium::P
    shape::S
end

CircleParticle{P, T, Dim} = Particle{T, Dim, P, Circle{T}}

volume(p::Particle) = volume(p.shape)

volume(particles::Vector{Particle}) = mapreduce(volume, +, particles)

function ==(p1::Particle, p2::Particle)
    p1.position == p2.position &&
    p1.medium == p2.medium &&
    p1.shape == p2.shape
end
