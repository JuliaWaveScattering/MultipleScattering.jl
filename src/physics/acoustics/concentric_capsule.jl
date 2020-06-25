
"""
    t_matrix(CapsuleParticle{T,2,Acoustic{T,2},Circle{T}}, Acoustic{T,2}, ω, order)

The T-matrix for a 2D circlular capsule particle in an acoustic medium.
"""
function t_matrix(cap::CapsuleParticle{T,2,Acoustic{T,2},Circle{T}}, medium::Acoustic{T,2}, ω::T, M::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    k = ω / medium.c
    k0 = ω / cap.inner.medium.c
    k1 = ω / cap.outer.medium.c
    a0 = outer_radius(cap.inner)
    a1 = outer_radius(cap.outer)
    q  = (medium.ρ*medium.c)/(cap.outer.medium.ρ*cap.outer.medium.c)
    q0 = (cap.inner.medium.ρ*cap.inner.medium.c)/(cap.outer.medium.ρ*cap.outer.medium.c)

    Yn(n::Integer)   = hankelh1(n,k1*a1)*besselj(n,k1*a0) - hankelh1(n,k1*a0)*besselj(n,k1*a1)
    Yddn(n::Integer) = diffhankelh1(n,k1*a1)*diffbesselj(n,k1*a0) - diffhankelh1(n,k1*a0)*diffbesselj(n,k1*a1)

    Ydn(n::Integer,x::Complex{T},y::Complex{T}) = hankelh1(n,x)*diffbesselj(n,y) - diffhankelh1(n,y)*besselj(n,x)

    function Tn(n::Integer)::Complex{T}
        numer = Ydn(n,k*a1,k*a1) * (q0*besselj(n,k0*a0)*Ydn(n,k1*a1,k1*a0) - Yn(n)*diffbesselj(n,k0*a0))
        denom = diffbesselj(n,k0*a0)*
            (q*hankelh1(n,k*a1)*Ydn(n,k1*a0,k1*a1) + diffhankelh1(n,k*a1)*Yn(n)) +
            q0*besselj(n,k0*a0)*
            (q*hankelh1(n,k*a1)*Yddn(n) - diffhankelh1(n,k*a1)*Ydn(n,k1*a1,k1*a0))

         return (numer - denom*besselj(n,k*a1)) / (denom*hankelh1(n,k*a1))
     end

    # Get Tns for positive m
    Tns = map(Tn,0:M)

    return Diagonal(vcat(reverse(Tns), Tns[2:end]))
end

function internal_field(x::AbstractArray{T}, p::CapsuleParticle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector) where T

    Nh = Int((length(scattering_coefficients) - one(T))/T(2.0)) #shorthand

    if iszero(p.outer.medium.c) || isinf(abs((p.outer.medium.c)))
        return zero(Complex{T})
    elseif norm(x - origin(p)) > outer_radius(p)
        @warn "point $x not insider particle $p. Returning zero."
        return zero(Complex{T})
    end

    k = ω / sim.source.medium.c
    k0 = ω / p.inner.medium.c
    k1 = ω / p.outer.medium.c
    a0 = outer_radius(p.inner)
    a1 = outer_radius(p.outer)
    q0 = (p.inner.medium.ρ*p.inner.medium.c) / (p.outer.medium.ρ*p.outer.medium.c)
    q = (sim.source.medium.ρ*sim.source.medium.c) / (p.outer.medium.ρ*p.outer.medium.c)

    Yn(n::Integer) = hankelh1(n,k1*a1)*besselj(n,k1*a0) - hankelh1(n,k1*a0)*besselj(n,k1*a1)
    Yddn(n::Integer) = diffhankelh1(n,k1*a1)*diffbesselj(n,k1*a0) - diffhankelh1(n,k1*a0)*diffbesselj(n,k1*a1)
    Ydn(n::Integer,x::Complex{T},y::Complex{T}) = hankelh1(n,x)*diffbesselj(n,y) - diffhankelh1(n,y)*besselj(n,x)

    denom(n::Integer) = q0 * besselj(n,a0*k0) *
            (q*besselj(n,a1*k)*Yddn(n) - diffbesselj(n,a1*k)*Ydn(n,a1*k1,a0*k1)) + diffbesselj(n,a0*k0) * (q*besselj(n,k*a1)*Ydn(n,a0*k1,a1*k1) + diffbesselj(n,k*a1)*Yn(n))
    forces = scattering_coefficients .* [Ydn(n,k*a1,k*a1)/denom(n) for n = -Nh:Nh];

    if norm(x - origin(p)) <= outer_radius(p.inner)
        coefs = - q0 * forces .* Ydn.(-Nh:Nh, a0*k1, a0*k1)
        basis = regular_basis_function(p.inner, ω)

        return sum(basis(Nh,x-origin(p)) .* coefs)
    else # calculate field in the mid region
        J_coefs = forces .* [
            q0*besselj(n, a0*k0)*diffhankelh1(n, a0*k1) - hankelh1(n, a0*k1)*diffbesselj(n, a0*k0)
        for n = -Nh:Nh]

        H_coefs = forces .* [
            besselj(n, a0*k1)*diffbesselj(n,a0*k0) - q0*besselj(n, a0*k0)*diffbesselj(n,a0*k1)
        for n = -Nh:Nh]

        J_basis = regular_basis_function(p.outer, ω)
        H_basis = outgoing_basis_function(p.outer.medium, ω)
        return sum(
            J_basis(Nh,x-origin(p)) .* J_coefs .+ H_basis(Nh,x-origin(p)) .* H_coefs
        )
    end
end
