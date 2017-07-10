
type Particle{T}
    x::Vector{T} #position of the centre
    r::T # radius
    c::Complex{T} # phase wave speed inside the particle
    ρ::T # density
end

function Particle{T}(x::Vector{T}; r::T = 1.0, c::Complex{T} = one(Complex{T}), ρ::T = zero(T))
    Particle{T}(x,r,c,ρ)
end

function Particle{T}(x::Vector{T}, r::T; c::Complex{T} = one(Complex{T}), ρ::T = zero(T))
    Particle{T}(x,r,c,ρ)
end

function volume{T}(particles::Vector{Particle{T}})
    mapreduce(volume, +, particles)
end

function volume{T}(p::Particle{T})
    π*(p.r)^2
end

"A kind of constructor for an array of homogenous particles"
function array_of_particles{T}(N::Int, a::T; c=one(Complex{T}), ρ=one(T))
    particles = Vector{Particle{T}}(N)
    for i=1:N
        particles[i] = Particle(zeros(T, 2), a, c, ρ)
    end
    return particles
end

function mean_radius{T}(particles::Vector{Particle{T}})
    radius_fnc(particle) = particle.r
    mapreduce(radius_fnc, +, particles) / length(particles)
end

function mean_volume{T}(particles::Vector{Particle{T}})
    volume(particles) / length(particles)
end

function isequal(p1::Particle, p2::Particle)
    (p1.x == p2.x) && (p1.r == p2.r) && (p1.c == p2.c) && (p1.ρ == p2.ρ)
end

==(p1::Particle, p2::Particle) = isequal(p1, p2)

!=(p1::Particle, p2::Particle) = !isequal(p1, p2)