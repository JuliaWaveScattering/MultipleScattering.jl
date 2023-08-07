"""
    t_matrix(Particle{2,Acoustic{T,2},Sphere{T,2}}, Acoustic{T,2}, ω, order)

The T-matrix for a 2D circlular Helmholtz resonator in a 2D acoustic medium.
"""
function t_matrix(p::Particle{2,Acoustic{T,2},SphericalHelmholtz{T,2}}, outer_medium::Acoustic{T,2}, ω::T, basis_order::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    M = basis_order
    ε = p.shape.aperture
    γe = 0.5772156649
    ak = outer_radius(p)*ω/outer_medium.c
    Qsum = sum([(besselj(m, ak)*diffhankelh1(m, ak) + diffbesselj(m, ak)*hankelh1(m, ak))^2/(diffbesselj(m, ak)*diffhankelh1(m, ak)) for m in -2M:2M])
    hε = (4im / pi) * (γe - 1im*pi/2 + log(ε/4)) - Qsum / 2

    # Check for material properties that don't make sense or haven't been implemented
    check_material(p, outer_medium)

    if p.medium.ρ < Inf || real(p.medium.c) < Inf
        @warn "Theory not done for general Helmholtz resonators, using Newmann boundaries instead."
    end

    "Returns a ratio used in multiple scattering which reflects the material properties of the particles"
    function Zn(m::Integer)::Complex{T}
        m = T(abs(m))

        # set the scattering strength and type
        numer = diffbesselj(m, ak)
        denom = diffhankelh1(m, ak)
        resonance_term = 2 / (hε * (pi * ak * diffhankelh1(m, ak))^2)

        return (numer / denom) - resonance_term
    end

    # Get Zns for positive m
    Zns = map(Zn,0:M)

    return - Diagonal(vcat(reverse(Zns), Zns[2:end]))
end
