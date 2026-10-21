
random_particles(particle_medium::PhysicalMedium{Dim}, particle_shape::Shape{Dim}; kws...) where {Dim} = random_particles(particle_medium, [particle_shape]; kws...)


function random_particles(particle_medium::PhysicalMedium{Dim}, particle_shapes::Vector{S};
        num_particles::Integer = length(particle_shapes), volume_fraction = 0,
        region_shape::Shape{Dim} = Box(fill(10*sum(outer_radius.(particle_shapes)),Dim)),
        current_particles::AbstractParticles{Dim} = AbstractParticle{Dim}[],
        kws...
) where {Dim, S<:Shape{Dim}}

    # Either the number of particles, or the volume fraction, for each shape
    amount_each = if volume_fraction == 0
        Int(round(num_particles / length(particle_shapes)))
    else
        volume_fraction / length(particle_shapes)
    end

    for s in particle_shapes
        add_particles = random_particles(particle_medium, s, region_shape, amount_each; current_particles = current_particles, kws...)
        current_particles = [current_particles; add_particles]
    end

    return current_particles
end

function random_particles(particle_medium::P, particle_shape::S,
    region_shape::Shape{Dim}, volfrac::AbstractFloat; kws...
) where {Dim,P<:PhysicalMedium{Dim},S<:Shape{Dim}}

    N = Int(round(volfrac * volume(region_shape) / volume(particle_shape)))
    return random_particles(particle_medium, particle_shape, region_shape, N; kws...)
end

"""
    random_particles(particle_medium, particle_shapes::Vector{Shape}, region_shape, volume_fraction::Number;
        seed=Random.make_seed())
    random_particles(particle_medium, particle_shape::Shape, region_shape, volume_fraction::Number;
        seed=Random.make_seed())
    random_particles(particle_medium, particle_shape::Shape, region_shape, N::Integer;
        seed=Random.make_seed())

Generate `N` random particles that fit inside `region_shape` (or fill with `volume_fraction`)

Specify seed to make output deterministic. Algorithm places particles unifomly randomly inside
the bounding box of `region_shape` and discards particle if it overlaps (based on outer radius)
or does not lies completely in box.

When passing particle_shapes::Vector{Shape} we assume each element is equally likely to occur. Repeating the same shape will lead to it being placed more often.
"""
function random_particles(particle_medium::P, particle_shape::S, region_shape::Shape{Dim}, N::Int;
        seed = -1,
        randgen = (seed == -1) ? Random.MersenneTwister() : Random.MersenneTwister(seed),
        verbose::Bool = false,
        separation_ratio = 1.0, # Min distance between particle centres relative to their outer radiuses.
        max_attempts_to_place_particle::Int = 3000, # Maximum number of attempts to place a particle
        current_particles::AbstractParticles{Dim} = AbstractParticle{Dim}[] # Particles already present.
) where {Dim,P<:PhysicalMedium{Dim},S<:Shape{Dim}}
    
    # Check volume fraction is not impossible
    volfrac = N * volume(particle_shape) / volume(region_shape)
    max_packing = if length(current_particles) > 0
        0.7854 - sum(volume.(current_particles))/volume(region_shape)
    else 0.7854
    end
    if volfrac  > max_packing
        error("Specified volume fraction is larger than optimal packing of circles.")
    end

    # Create pseudorandom device with specified seed

    box = bounding_box(region_shape)

    if verbose
        @printf("""\n
            Generating %d randomly positioned %s shaped particles
            Total particle volume: %0.5g
            Inside %s of volume: %0.5g
            Particle volume fraction: %0.5g
            Bounding box volume: %0.5g
            """, N, name(particle_shape), N*volume(particle_shape), name(region_shape),
            volume(region_shape), volfrac, volume(box)
        )
    end

    # Allocate memory for the new particles. These can have a different type from current_particles.
    particles = Vector{Particle{Dim,P,S}}(undef, N)

    too_close(p1, p2) = norm(origin(p1) - origin(p2)) < separation_ratio * (outer_radius(p1) + outer_radius(p2))

    for n = 1:N

        num_attempts = 0
        overlapping = true
        while overlapping

            outside_box = true
            while outside_box
                x = (box.dimensions ./ 2) .* (1 .- 2 .* rand(randgen,typeof(origin(box)))) + origin(box)
                if S <: Helmholtz
                    θ = 2π*rand(randgen)
                    particles[n] = Particle(particle_medium, congruent(particle_shape, x, θ))
                else
                    particles[n] = Particle(particle_medium, congruent(particle_shape, x))
                end
                outside_box = !(particles[n] ⊆ region_shape)
            end

            # compare with the current particles and the previously placed particles
            overlapping = any(p -> too_close(particles[n], p), current_particles) ||
                any(i -> too_close(particles[n], particles[i]), 1:(n-1))

            num_attempts += 1
            if num_attempts > max_attempts_to_place_particle
                error("Tried to place a scatterer $max_attempts_to_place_particle times unsuccessfully. You could try increasing max_attempts_to_place_particle")
            end

        end
    end

    return particles
end
