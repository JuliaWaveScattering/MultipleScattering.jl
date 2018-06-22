
type Particle{T}
    x::Vector{T} #position of the centre
    r::T # radius
    c::Complex{T} # phase wave speed inside the particle
    ρ::T # the particles density
end

function Particle{T}(x::Vector{T}, r::T = 1.0; c::Complex{T} = one(Complex{T}), ρ::T = zero(T))
    Particle{T}(x,r,c,ρ)
end

function volume{T}(particles::Vector{Particle{T}})
    isempty(particles) ? zero(T) : mapreduce(volume, +, particles)
end

function volume{T}(p::Particle{T})
    π*(p.r)^2
end

function mean_radius{T}(particles::Vector{Particle{T}})
    radius_fnc(particle) = particle.r
    isempty(particles) ? zero(T) : mapreduce(radius_fnc, +, particles) / length(particles)
end

function std_radius{T}(particles::Vector{Particle{T}})
    radii = [p.r for p in particles]
    std(radii)
end

function isequal(p1::Particle, p2::Particle)
    isequal(p1.x,p2.x) && isequal(p1.r,p2.r) &&
    isequal(p1.c,p2.c) && isequal(p1.ρ,p2.ρ)
end

function ==(p1::Particle, p2::Particle)
    (p1.x == p2.x) && (p1.r == p2.r) && (p1.c == p2.c) && (p1.ρ == p2.ρ)
end
