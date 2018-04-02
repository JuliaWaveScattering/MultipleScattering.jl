
abstract PhysicalProperties{T,Dim,FieldDim} where T <: AbstractFloat, Dim::Integer, FieldDim::Integer end

type Particle{T,Dim,P,S} where P <: PhysicalProperties{T,Dim,FieldDim}, S <: Shape{T, Dim}
    position::MVector{T,Dim}
    medium::P
    shape::S
end

CircleParticle{P, T, Dim} = Particle{T, Dim, P, Circle{T, Dim}}
