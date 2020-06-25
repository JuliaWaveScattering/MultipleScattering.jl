AcousticCircleParticle{T} = Particle{T,2,Acoustic{T,2},Circle{T}}

"""
    t_matrix(Particle{T,2,Acoustic{T,2},Circle{T}}, Acoustic{T,2}, ω, order)

The T-matrix for a 2D circlular acoustic particle in a 2D acoustic medium.
"""
function t_matrix(p::Particle{T,2,Acoustic{T,2},Circle{T}}, outer_medium::Acoustic{T,2}, ω::T, basis_order::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    # Check for material properties that don't make sense or haven't been implemented
    check_material(p, outer_medium)

    "Returns a ratio used in multiple scattering which reflects the material properties of the particles"
    function Zn(m::Integer)::Complex{T}
        m = T(abs(m))
        ak = outer_radius(p)*ω/outer_medium.c

        # set the scattering strength and type
        if isinf(p.medium.c) || isinf(p.medium.ρ)
            numer = diffbesselj(m, ak)
            denom = diffhankelh1(m, ak)
        elseif iszero(outer_medium.ρ)
            numer = diffbesselj(m, ak)
            denom = diffhankelh1(m, ak)
        elseif iszero(p.medium.ρ) || iszero(p.medium.c)
            numer = besselj(m, ak)
            denom = hankelh1(m, ak)
        else
            q = impedance(p.medium)/impedance(outer_medium) # Impedance ratio
            γ = outer_medium.c/p.medium.c #speed ratio
            numer = q * diffbesselj(m, ak) * besselj(m, γ * ak) - besselj(m, ak)*diffbesselj(m, γ * ak)
            denom = q * diffhankelh1(m, ak) * besselj(m, γ * ak) - hankelh1(m, ak)*diffbesselj(m, γ * ak)
        end

        return numer / denom
    end

    # Get Zns for positive m
    Zns = map(Zn,0:basis_order)

    return - Diagonal(vcat(reverse(Zns), Zns[2:end]))
end

"""
    internal_field(x::SVector{2,T}, p::Particle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}})

The internal field for a 2D circlular acoustic particle in a 2D acoustic medium.
"""
function internal_field(x::AbstractVector{T}, p::Particle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}}) where T

    Nh = Int((length(scattering_coefficients) - one(T))/T(2.0)) #shorthand
    if iszero(p.medium.c) || isinf(abs(p.medium.c))
        return zero(Complex{T})
    else
        r = outer_radius(p)
        k = ω/sim.source.medium.c
        kp = ω/p.medium.c
        diagZ = - diag(t_matrix(p, sim.source.medium, ω, Nh))

        internal_coefs = scattering_coefficients ./
            (diagZ .* besselj.(-Nh:Nh,kp*r)) .*
            (diagZ .* hankelh1.(-Nh:Nh,k*r) - besselj.(-Nh:Nh,k*r))

        inner_basis = regular_basis_function(p, ω)

        return sum(inner_basis(Nh, x-origin(p)) .* internal_coefs)
    end
end
