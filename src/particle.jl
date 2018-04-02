
type Particle{T,Dim,P,S} where P <: PhysicalProperties{Dim,FieldDim,T}, S <: Shape{Dim, T}
    position::MVector{Dim,T}
    medium::P
    shape::S
end

CircleParticle{P, T, Dim} = Particle{T, Dim, P, Circle{T}}
