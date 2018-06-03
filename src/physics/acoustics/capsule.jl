
struct AcousticCapsule{T,Dim} <: PhysicalProperties{T,Dim,1}
    in::Acoustic{T,Dim} # inner properties
    out::Acoustic{T,Dim} # outer properties
    shape_in::T # inner homogenous shape
    shape_out::T # outer homogenous shape
end

name(a::AcousticCapsule{T,Dim}) where {Dim,T} = "$(Dim)D Acoustic Capsule"
