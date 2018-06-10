AcousticCircleParticle{T} = Particle{T,2,Acoustic{T,2},Circle{T}}

# T-matrix for a 2D circlular acoustic particle in a 2D acoustic medium
function t_matrix(p::Particle{T,2,Acoustic{T,2},Circle{T}}, outer_medium::Acoustic{T,2}, ω::T, M::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    # Check for material properties that don't make sense or haven't been implemented
    if isnan(abs(p.medium.c)*p.medium.ρ)
        throw(DomainError("Particle's phase speed times density is not a number!"))
    elseif isnan(abs(outer_medium.c)*outer_medium.ρ)
        throw(DomainError("The medium's phase speed times density is not a number!"))
    elseif iszero(outer_medium.c)
        throw(DomainError("Wave propagation in a medium with zero phase speed is not defined"))
    elseif iszero(outer_medium.ρ) && iszero(p.medium.c*p.medium.ρ)
        throw(DomainError("Scattering in a medium with zero density from a particle with zero density or zero phase speed is not defined"))
    elseif iszero(outer_radius(p))
        throw(DomainError("Scattering from a circle of zero radius is not implemented yet"))
    end

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
            q = (p.medium.c*p.medium.ρ)/(outer_medium.c*outer_medium.ρ) #the impedance
            γ = outer_medium.c/p.medium.c #speed ratio
            numer = q * diffbesselj(m, ak) * besselj(m, γ * ak) - besselj(m, ak)*diffbesselj(m, γ * ak)
            denom = q * diffhankelh1(m, ak) * besselj(m, γ * ak) - hankelh1(m, ak)*diffbesselj(m, γ * ak)
        end

        return numer / denom
    end

    # Get Zns for positive m
    Zns = map(Zn,0:M)

    return - Diagonal(vcat(reverse(Zns), Zns[2:end]))
end


function internal_field(x::SVector{2,T}, p::Particle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}}) where T

    Nh = Int((length(scattering_coefficients) - one(T))/T(2.0)) #shorthand
    if iszero(p.medium.c) || isinf(abs(p.medium.c))
        return zero(Complex{T})
    else
        r = outer_radius(p)
        k = ω/sim.medium.c
        kp = ω/p.medium.c
        Z = - t_matrix(p, sim.medium, ω, Nh)
        internal_coef(m::Int) = scattering_coefficients[m+Nh+1] / (Z[m+Nh+1,m+Nh+1]*besselj(m,kp*r)) * (Z[m+Nh+1,m+Nh+1]*hankelh1(m,k*r) - besselj(m,k*r))

        inner_basis = basis_function(p, ω)
        return sum(-Nh:Nh) do m
            inner_basis(m, x-origin(p)) * internal_coef(m)
        end
    end
end
