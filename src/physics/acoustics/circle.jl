AcousticCircleParticle{T} = Particle{T,2,Acoustic{T,2},Circle{T}}

function inner_basis_coefficients(p::Particle{T,2,Acoustic{T,2},Circle{T}}, medium::Acoustic{T,2}, ω::T, scattering_coefficients::AbstractVector; basis_order::Int=5) where T

    Nh = basis_order
    if iszero(p.medium.c) || isinf(abs(p.medium.c))
        return zeros(Complex{Float64},2Nh+1)
    else
        r = outer_radius(p)
        k = ω/medium.c
        kp = ω/p.medium.c
        Z = - t_matrix(p.shape, p.medium, medium, ω, basis_order)
        return map(-Nh:Nh) do m
             scattering_coefficients[m+Nh+1] / (Z[m+Nh+1,m+Nh+1]*besselj(m,kp*r)) * (Z[m+Nh+1,m+Nh+1]*hankelh1(m,k*r) - besselj(m,k*r))
        end
    end
end

# T-matrix for a 2D circlular acoustic particle in a 2D acoustic medium
function t_matrix(circle::Circle{T}, inner_medium::Acoustic{T,2}, outer_medium::Acoustic{T,2}, ω::T, M::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    # Check for material properties that don't make sense or haven't been implemented
    if isnan(abs(inner_medium.c)*inner_medium.ρ)
        throw(DomainError("Particle's phase speed times density is not a number!"))
    elseif isnan(abs(outer_medium.c)*outer_medium.ρ)
        throw(DomainError("The medium's phase speed times density is not a number!"))
    elseif iszero(outer_medium.c)
        throw(DomainError("Wave propagation in a medium with zero phase speed is not defined"))
    elseif iszero(outer_medium.ρ) && iszero(inner_medium.c*inner_medium.ρ)
        throw(DomainError("Scattering in a medium with zero density from a particle with zero density or zero phase speed is not defined"))
    elseif iszero(circle.radius)
        throw(DomainError("Scattering from a circle of zero radius is not implemented yet"))
    end

    "Returns a ratio used in multiple scattering which reflects the material properties of the particles"
    function Zn(m::Integer)::Complex{T}
        m = T(abs(m))
        ak = circle.radius*ω/outer_medium.c

        # set the scattering strength and type
        if isinf(inner_medium.c) || isinf(inner_medium.ρ)
            numer = diffbesselj(m, ak)
            denom = diffhankelh1(m, ak)
        elseif iszero(outer_medium.ρ)
            numer = diffbesselj(m, ak)
            denom = diffhankelh1(m, ak)
        elseif iszero(inner_medium.ρ) || iszero(inner_medium.c)
            numer = besselj(m, ak)
            denom = hankelh1(m, ak)
        else
            q = (inner_medium.c*inner_medium.ρ)/(outer_medium.c*outer_medium.ρ) #the impedance
            γ = outer_medium.c/inner_medium.c #speed ratio
            numer = q * diffbesselj(m, ak) * besselj(m, γ * ak) - besselj(m, ak)*diffbesselj(m, γ * ak)
            denom = q * diffhankelh1(m, ak) * besselj(m, γ * ak) - hankelh1(m, ak)*diffbesselj(m, γ * ak)
        end

        return numer / denom
    end

    # Get Zns for positive m
    Zns = map(Zn,0:M)

    return - Diagonal(vcat(reverse(Zns), Zns[2:end]))
end
