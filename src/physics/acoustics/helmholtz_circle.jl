"""
    t_matrix(Particle{2,Acoustic{T,2},Sphere{T,2}}, Acoustic{T,2}, ω, order)

The T-matrix for a 2D circlular Helmholtz resonator in a 2D acoustic medium.
"""
function t_matrix(p::Particle{2,Acoustic{T,2},SphericalHelmholtz{T,2}}, outer_medium::Acoustic{T,2}, ω::T, basis_order::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    M = basis_order
    ε = p.shape.aperture
    γe = 0.5772156649
    k = ω/outer_medium.c
    kb = k*outer_radius(p)
    Qsum = sum([(besselj(m, kb)*diffhankelh1(m, kb) + diffbesselj(m, kb)*hankelh1(m, kb))^2/(diffbesselj(m, kb)*diffhankelh1(m, kb)) for m in -2M:2M])
    hε = (4im / pi) * (γe - 1im*pi/2 + log(k*ε/4)) - Qsum / 2

    # Check for material properties that don't mkbe sense or haven't been implemented
    check_material(p, outer_medium)

    if p.medium.ρ < Inf || real(p.medium.c) < Inf
        @warn "Theory not done for general Helmholtz resonators, using Newmann boundaries instead."
    end

    "Returns a ratio used in multiple scattering which reflects the material properties of the particles"
    function Zn(m::Integer)::Complex{T}
        m = T(abs(m))

        # set the scattering strength and type
        numer = diffbesselj(m, kb)
        denom = diffhankelh1(m, kb)
        resonance_term = 2 / (hε * (pi * kb * diffhankelh1(m, kb))^2)

        return (numer / denom) - resonance_term
    end

    # Get Zns for positive m
    Zns = map(Zn,0:M)

    return - Diagonal(vcat(reverse(Zns), Zns[2:end]))
end
