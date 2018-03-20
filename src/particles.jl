
type Particles{P,T,Dim} where P <: PhysicalProperties{T}
    position::MVector{T,Dim}
    medium::P
    shape::Shape{T}
end
