# T-matrix for a 2D circlular acoustic particle in a 2D acoustic medium

function inner_basis_coefficients(x::SVector{2,T}, p::CapsuleParticle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector;
        basis_order::Int=5) where T

    medium = sim.medium
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

function t_matrix(cap::CapsuleParticle{T,2,Acoustic{T,2},Circle{T}}, medium::Acoustic{T,2}, ω::T, M::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    k = ω / medium.c
    k0 = ω / cap.inner.medium.c
    k1 = ω / cap.outer.medium.c
    a0 = outer_radius(cap.inner)
    a1 = outer_radius(cap.outer)
    q  = (medium.ρ*medium.c)/(cap.outer.medium.ρ*cap.outer.medium.c)
    q0 = (cap.inner.medium.ρ*cap.inner.medium.c)/(cap.outer.medium.ρ*cap.outer.medium.c)

    Yn(n::Integer) = hankelh1(n,k1*a1)*besselj(n,k1*a0) - hankelh1(n,k1*a0)*besselj(n,k1*a1)
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

    return - Diagonal(vcat(reverse(Tns), Tns[2:end]))
end
