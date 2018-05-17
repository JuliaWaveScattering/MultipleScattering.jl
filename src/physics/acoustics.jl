
"""
Physical properties for a homogenous isotropic acoustic medium. Produces a
scalar (1D) field in arbitrary dimensions.
"""
type Acoustic{Dim,T} <: PhysicalProperties{Dim,1,T}
    ρ::T # Density
    c::Complex{T} # Phase velocity
end

# Constructor which supplies the dimension without explicitly mentioning type
Acoustic(ρ::T,c::Complex{T},Dim::Integer) where {T} =  Acoustic{Dim,T}(ρ,c)

name(a::Acoustic{Dim,T}) where {Dim,T} = "$(Dim)D Acoustic"

function basis_function(medium::Acoustic{2,T}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::SVector{2,T})
        r = norm(x)
        θ = atan2(x[2],x[1])
        k = ω/medium.c
        hankelh1(m,k*r)*exp(im*θ*m)
    end
end

"Basis function when inside a particle. Assumes particle is a circle, which approximately works for all shapes."
function basis_function(p::Particle{2,Acoustic{2,T}}, ω::T) where {T}
    return function acoustic_basis_function(m::Integer, x::SVector{2,T})
        r = norm(x)
        θ = atan2(x[2],x[1])
        k = ω/p.medium.c
        besselj(m,k*r)*exp(im*θ*m)
    end
end



# Type aliases for convenience
TwoDimAcoustic{T} = Acoustic{2,T}
AcousticCircleParticle{T} = Particle{2, Acoustic{2, T}, Circle{T}, T}
TwoDimAcousticFrequencySimulation{T} = FrequencySimulation{2,Acoustic{2,T},T}


"""
Physical properties for a capsule with the same inner and outer shape, where the
characteristic length scale of the two shapes is `size_ratio`. Inner and outer
are both a homogenous isotropic acoustic medium. Produces a scalar (1D) field in
arbitrary dimensions.
"""
type AcousticCapsule{T,Dim} <: PhysicalProperties{Dim,1,T}
    inner_ρ::T # Density in inner shape
    inner_c::Complex{T} # Phase velocity in inner shape
    outer_ρ::T # Density in outer shape
    outer_c::Complex{T} # Phase velocity in outer shape
    size_ratio::T # Ratio of size of inner shape to outer shape
end

name(a::AcousticCapsule{T,Dim}) where {Dim,T} = "$(Dim)D Acoustic Capsule"


# T-matrix for a 2D circlular acoustic particle in a 2D acoustic medium
function t_matrix(circle::Circle{T}, inner_medium::Acoustic{2,T}, outer_medium::Acoustic{2,T}, ω::T, M::Integer)::Diagonal{Complex{T}} where T<:AbstractFloat

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

    return - Diagonal(vcat(reverse(Zns), Zns[2:end]))
end


function TwoDimAcousticPointSource{T}(medium::Acoustic{2,T}, source_position::AbstractVector{T}, amplitude::T = one(T))::Source{Acoustic{2,T},T}
    field(x,ω) = amplitude*Complex{T}(im/4)*hankelh1(0,ω/medium.c*norm(x-source_position))
    # using Graf's addition theorem
    coef(n,centre,ω) = amplitude*Complex{T}(im/4)*hankelh1(-n,ω/medium.c*norm(centre - source_position))*
exp(-Complex{T}(im)*n*atan2(centre[2] - source_position[2], centre[1] - source_position[1]))

    return Source{Acoustic{2,T},T}(field,coef)
end

function TwoDimAcousticPlanarSource{T}(medium::Acoustic{2,T}, source_position::AbstractVector{T}, source_direction::AbstractVector{T} = [one(T),zero(T)], amplitude::T = one(T))::Source{Acoustic{2,T},T}
    field(x,ω) = amplitude*exp(im*ω/medium.c*dot(x-source_position,source_direction))
    # Jacobi-Anger expansion
    coef(n,centre,ω) = field(centre,ω) * exp(im * n *(T(pi)/2 - atan2(source_direction[2],source_direction[1]) ))

    return Source{Acoustic{2,T},T}(field,coef)
end

function inner_basis_coefficients(p::Particle{2, Acoustic{2,T}}, medium::Acoustic{2,T}, ω::T, scattering_coefficients::AbstractVector; basis_order::Int=5) where T<:Number
    r = outer_radius(p)
    k = ω/medium.c
    kp = ω/p.medium.c
    Nh = basis_order

    Z = - t_matrix(p.shape, p.medium, medium, ω, basis_order)
    return map(-Nh:Nh) do m
         scattering_coefficients[m+Nh+1] / (Z[m+Nh+1,m+Nh+1]*besselj(m,kp*r)) * (Z[m+Nh+1,m+Nh+1]*hankelh1(m,k*r) - besselj(m,k*r))
    end
end

function besselj_field(source::Source{Acoustic{2,T},T}, medium::Acoustic{2,T}, centre::AbstractVector{T}; basis_order = 4) where T<:Number

    field(x,ω) = sum(
        source.coef(n,centre,ω)*besselj(n,ω/medium.c*norm(x - centre))*exp(im*n*atan2(x[2] - centre[2],x[1] - centre[1]))
    for n = -basis_order:basis_order)

   return field
end
