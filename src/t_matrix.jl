

function t_matrix(circle::Circle{T}, inner_medium::Acoustic{2,T}, outer_medium::Acoustic{2,T}, ω::T, M::Integer)::Diagonal{T} where T<:AbstractFloat

    # Check for material properties that don't make sense or haven't been implemented
    if isnan(inner_medium.c*inner_medium.ρ)
        throw(DomainError("Scattering from a particle with zero density or zero phase speed is not defined"))
    elseif isnan(outer_medium.c*outer_medium.ρ)
        throw(DomainError("Wave propagation in a medium with zero density or zero phase speed is not defined"))
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
            γ = outer_medium.c/inner_medium.c #speed ratio
            numer = diffbesselj(m, ak) * besselj(m, γ * ak)
            denom = diffhankelh1(m, ak) * besselj(m, γ * ak)
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

    return Diagonal(vcat(reverse(Zns), Zns[2:end]))
end

"""
Returns a 2M+1 by 2M+1 T-matrix for particle with specific shape, physical
properties in a medium with a specific physical property at a specific angular
wavenumber.
"""
function t_matrix(shape::Shape{T}, inner_medium::PhysicalProperties{Dim,FieldDim,T}, outer_medium::PhysicalProperties{Dim,FieldDim,T}, ω::T, M::Integer)::AbstractMatrix{T} where {Dim,FieldDim,T<:AbstractFloat}
    error("T-matrix function is not yet written for $(name(inner_medium)) $(name(shape)) in a $(name(outer_medium)) medium")
end
