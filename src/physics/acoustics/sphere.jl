"""
    t_matrix(Particle{T,3,Acoustic{T,3},Sphere{T}}, Acoustic{T,3}, ω, order)

The T-matrix for a 3D spherical acoustic particle in a 3D acoustic medium.
"""
function t_matrix(p::Particle{T,3,Acoustic{T,3},Sphere{T}}, outer_medium::Acoustic{T,3}, ω::T, basis_order::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

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

    # Get Zns for positive m
    Zns = map(Zn,0:basis_order)

    return - Diagonal(Zns)
end

# """
#     internal_field(x::SVector{2,T}, p::Particle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}})
#
# The internal field for a 2D circlular acoustic particle in a 2D acoustic medium.
# """
# function internal_field(x::SVector{2,T}, p::Particle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}}) where T

    # Nh = Int((length(scattering_coefficients) - one(T))/T(2.0)) #shorthand
    # if iszero(p.medium.c) || isinf(abs(p.medium.c))
        # return zero(Complex{T})
    # else
    #     r = outer_radius(p)
    #     k = ω/sim.source.medium.c
    #     kp = ω/p.medium.c
    #     diagZ = - diag(t_matrix(p, sim.source.medium, ω, Nh))
    #
    #     internal_coefs = scattering_coefficients ./
    #         (diagZ .* besselj.(-Nh:Nh,kp*r)) .*
    #         (diagZ .* hankelh1.(-Nh:Nh,k*r) - besselj.(-Nh:Nh,k*r))
    #
    #     inner_basis = regular_basis_function(p, ω)
    #
    #     return sum(inner_basis(Nh, x-origin(p)) .* internal_coefs)
    # end
# end
