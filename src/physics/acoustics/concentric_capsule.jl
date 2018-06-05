# T-matrix for a 2D circlular acoustic particle in a 2D acoustic medium

function internal_field(x::SVector{2,T}, p::CapsuleParticle{T,2,Acoustic{T,2},Circle{T}}, sim::FrequencySimulation{T,2,Acoustic{T,2}}, ω::T, scattering_coefficients::AbstractVector;
        basis_order::Int=5) where T

    Nh = basis_order

    if iszero(p.outer.medium.c) || isinf(abs((p.outer.medium.c)))
        return zero(Complex{T})
    elseif norm(x - origin(p)) > outer_radius(p)
        warn("point $x not insider particle $p. Returning zero.")
        return zero(Complex{T})
    end

    k = ω / medium.c
    k0 = ω / p.inner.medium.c
    k1 = ω / p.outer.medium.c
    a0 = outer_radius(p.inner)
    a1 = outer_radius(p.outer)
    q0 = (p.inner.medium.ρ*p.inner.medium.c)/(p.outer.medium.ρ*p.outer.medium.c)

    Yn(n::Integer) = hankelh1(n,k1*a1)*besselj(n,k1*a0) - hankelh1(n,k1*a0)*besselj(n,k1*a1)
    Ydn(n::Integer,x::Complex{T},y::Complex{T}) = hankelh1(n,x)*diffbesselj(n,y) - diffhankelh1(n,y)*besselj(n,x)

    force(n) = scattering_coefficients[m+Nh+1] * hankelh1(n, k*a1) +
        sim.source.coef(n,origin(p),ω) * besselj(n, k*a1);

    if norm(x - origin(p)) <= outer_radius(p.inner)
        function coef(n::Int)
            numer = q0*Ydn(n, a0*k1, a0*k1)
            denom = q0*besselj(n, a0*k0)*Ydn(n, a1*k1, a0*k1) - diffbesselj(n, a0*k0)*Yn(n)
            return force(n) * numer / denom
        end
        basis = basis_function(p.inner, ω)
        return map(-Nh:Nh) do m
             basis(m,x) * coef(n)
        end
    end
    # else calculate field in mid region
    function J_coef(n::Int)
        numer = q0*besselj(n, a0*k0)*diffhankelh1(n, a0*k1) - hankelh1(n, a0*k1)*diffbesselj(n, a0*k0)
        denom = Yn(n)*diffbesselj(n, a0*k0) - q0*besselj(n, a0*k0)*Ydn(n, a1*k1, a0*k1)
        return force(n) * numer / denom
    end
    function H_coef(n::Int)
        numer =  besselj(n, a0*k1)*diffbesselj(n,a0*k0) - q0*besselj(n, a0*k0)*diffbesselj(n,a0*k1)
        denom = Yn(n)*diffbesselj(n,a0*k0) - q0*besselj(n, a0*k0)*Ydn(n, a1*k1, a0*k1)
        return force(n) * numer / denom
    end

    J_basis = basis_function(p.outer, ω)
    H_basis = basis_function(p.outer.medium, ω)
    # return map(-Nh:Nh) do m
    #      basis(m,x) * coef(n)
    # end

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
