"""
    t_matrix(Particle{3,Acoustic{T,3},Sphere{T,3}}, Acoustic{T,3}, ω, order)

The T-matrix for a 3D spherical acoustic particle in a 3D acoustic medium.
"""
function t_matrix(p::Particle{3,Acoustic{T,3},Sphere{T,3}}, outer_medium::Acoustic{T,3}, ω::T, basis_order::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    # Check for material properties that don't make sense or haven't been implemented
    check_material(p, outer_medium)

    "Returns a ratio used in multiple scattering which reflects the material properties of the particles"
    function Zn(m::Integer)::Complex{T}
        m = T(abs(m))
        ak = outer_radius(p)*ω/outer_medium.c

        # set the scattering strength and type
        if isinf(p.medium.c) || isinf(p.medium.ρ) || iszero(outer_medium.ρ)
            numer = diffsbesselj(m, ak)
            denom = diffshankelh1(m, ak)
        elseif iszero(p.medium.ρ) || iszero(p.medium.c)
            numer = sbesselj(m, ak)
            denom = shankelh1(m, ak)
        else
            q = impedance(p.medium)/impedance(outer_medium) # Impedance ratio
            γ = outer_medium.c/p.medium.c #speed ratio
            numer = q * diffsbesselj(m, ak) * sbesselj(m, γ * ak) - sbesselj(m, ak) * diffsbesselj(m, γ * ak)
            denom = q * diffshankelh1(m, ak) * sbesselj(m, γ * ak) - shankelh1(m, ak) * diffsbesselj(m, γ * ak)
        end

        return numer / denom
    end
    Zns = Zn.(0:basis_order)

    # the T-matrix is the same for every m with the same degree l
    T_vec = - [Zns[l+1] for l = 0:basis_order for m = -l:l]

    return Diagonal(T_vec)
end
