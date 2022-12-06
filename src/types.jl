"""
    AbstractSymmetry

An abstract types which dictates the symmetry of the setup. That is, the symmetry shared between the incident wave and the shape of the material.
"""
abstract type AbstractSymmetry{Dim} end

"""
    WithoutSymmetry

A type used to describe materials and incident waves which share no common symmetry. This will lead to the most general regular wave expansion for the eignvectors.
"""
struct WithoutSymmetry{Dim} <: AbstractSymmetry{Dim} end

"""
An incident plane-wave and infinite cylinder material will result in all fields being cylindrical Bessel waves with a phase.
"""
abstract type AbstractTranslationSymmetry{Dim} <: AbstractSymmetry{Dim} end

"""
    TranslationSymmetry

A type used to describe materials and incident waves that both share a translation symmetry.
"""
struct TranslationSymmetry{Dim,T} <: AbstractTranslationSymmetry{Dim} 
    direction::SVector{Dim,T}
end

"""
An incident plane-wave and halfspace material will result in all fields being plane-waves.
"""
abstract type AbstractPlanarSymmetry{Dim} <: AbstractSymmetry{Dim} end

"""
    PlanarSymmetry

A type used to describe materials and incident waves that both share a planar symmetry.
"""
struct PlanarSymmetry{Dim} <: AbstractPlanarSymmetry{Dim} end

"""
For spatial dimension > 2, we can consider problems that have azimuthal symmetry. For example, a plane-wave incident on a sphere.
"""
abstract type AbstractAzimuthalSymmetry{Dim} <: AbstractSymmetry{Dim} end

"""
    AzimuthalSymmetry

A type used to describe materials and incident waves in 3 spatial dimensions that share a symmetry about the azimuth.
"""
struct AzimuthalSymmetry{Dim} <: AbstractAzimuthalSymmetry{Dim} end
AzimuthalSymmetry() = AzimuthalSymmetry{3}()

"""
    RadialSymmetry

A type used to describe materials and incident waves in that are both radially symmetric. That is, the material is a sphere, and the incident wave is radially symmetric around the center of this spherical material.
"""
struct RadialSymmetry{Dim} <: AbstractAzimuthalSymmetry{Dim} end

"""
    PlanarAzimuthalSymmetry

A type used to describe a materail and incident wave which have both [`PlanarSymmetry`](@ref) and [`AzimuthalSymmetry`](@ref).
"""
struct PlanarAzimuthalSymmetry{Dim} <: AbstractPlanarSymmetry{Dim} end
PlanarAzimuthalSymmetry() = PlanarAzimuthalSymmetry{3}()

function azimuthalnormal(dim::Int)
    if dim == 2
        return [1,0]
    elseif dim == 3
        return [0,0,1]
    else
        return [zeros(dim-1);1]
    end
end

# If T1 is a shape, source, or material, first calculate its symmetry
"""
    Symmetry(T)

returns the symmetry of the object T. Some possible symmetries include [`PlanarSymmetry`](@ref),  [`RadialSymmetry`](@ref),  and [`WithoutSymmetry`](@ref).
"""


"""
    Symmetry(T1,T2)

returns the combined symmetry of the two objects T1 and T2.
"""
Symmetry(T1,T2) = Symmetry(Symmetry(T1),Symmetry(T2))

Symmetry(sym1::AbstractSymmetry{Dim},sym2::AbstractSymmetry{Dim}) where Dim = WithoutSymmetry{Dim}()
Symmetry(sym1::S,sym2::S) where S<:AbstractSymmetry = sym1

function Symmetry(sym1::TranslationSymmetry{Dim},sym2::TranslationSymmetry{Dim}) where Dim
    if abs(dot(sym1.direction,sym2.direction)) > 0
        PlanarSymmetry{Dim}()
    else sym1
    end    
end    

Symmetry(sym1::AbstractPlanarSymmetry{Dim},sym2::AbstractPlanarSymmetry{Dim}) where Dim = PlanarSymmetry{Dim}()
Symmetry(sym1::PlanarAzimuthalSymmetry{Dim},sym2::PlanarAzimuthalSymmetry{Dim}) where Dim = PlanarAzimuthalSymmetry{Dim}()

Symmetry(sym1::AzimuthalSymmetry{Dim},sym2::PlanarAzimuthalSymmetry{Dim}) where Dim = AzimuthalSymmetry{Dim}()
Symmetry(sym1::PlanarAzimuthalSymmetry{Dim},sym2::AzimuthalSymmetry{Dim}) where Dim = AzimuthalSymmetry{Dim}()

Symmetry(sym1::AbstractAzimuthalSymmetry{Dim},sym2::RadialSymmetry{Dim}) where Dim = AzimuthalSymmetry{Dim}()
Symmetry(sym1::RadialSymmetry{Dim},sym2::AbstractAzimuthalSymmetry{Dim}) where Dim = AzimuthalSymmetry{Dim}()

Symmetry(sym1::RadialSymmetry{Dim},sym2::RadialSymmetry{Dim}) where Dim = RadialSymmetry{Dim}()

Symmetry(sym1::PlanarAzimuthalSymmetry{Dim},sym2::RadialSymmetry{Dim}) where Dim = AzimuthalSymmetry{Dim}()
Symmetry(sym1::RadialSymmetry{Dim},sym2::PlanarAzimuthalSymmetry{Dim}) where Dim = AzimuthalSymmetry{Dim}()
