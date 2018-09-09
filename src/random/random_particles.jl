const MAX_ATTEMPTS_TO_FIT_PARTICLE = 3000

"""
    random_particles(particle_medium, particle_shape, box_shape, volume_fraction::Number;
        seed=Base.Random.make_seed())
    random_particles(particle_medium, particle_shape, box_shape, N::Integer;
        seed=Base.Random.make_seed())

Generate `N` random particles that fit inside `box_shape` (or fill with `volume_fraction`)

Specify seed to make output deterministic. Algorithm places particles unifomly randomly inside
the bounding rectangle of `box_shape` and discards particle if it overlaps (based on outer radius)
or does not lies completely in box box.
"""
function random_particles(particle_medium::P, particle_shape::S,
    box_shape::Shape{T,Dim}, N::Integer; seed=Base.Random.make_seed()) where {T,Dim,P<:PhysicalProperties{T,Dim},S<:Shape{T,Dim}}

    # Check volume fraction is not impossible
    volfrac = N * volume(particle_shape) / volume(box_shape)
    if volfrac > 0.7854
        error("Specified volume fraction is larger than optimal packing of circles.")
    end

    # Create pseudorandom device with specified seed
    randgen = MersenneTwister(seed)

    # Min distance between the centre of particles relative to their outer radiuses
    # The factor 2.1radius is seen in several papers (therefore ratio of 1.05)
    separation_ratio = 1.05

    bounding_rect = bounding_rectangle(box_shape)
    bounding_rect_size = SVector(2*bounding_rect.width, 2*bounding_rect.height)

    @printf("""\n
        Generating %d randomly positioned %s shaped particles
        Total particle volume: %0.5g
        Inside %s of volume: %0.5g
        Particle volume fraction: %0.5g
        Bounding box volume: %0.5g
        """, N, name(particle_shape), N*volume(particle_shape), name(box_shape),
        volume(box_shape), volfrac, volume(bounding_rect)
    )

    # Allocate memory for particles
    particles = Vector{Particle{T,Dim,P,S}}(N)

    for n = 1:N

        num_attempts = 0
        overlapping = true
        while overlapping

            outside_box = true
            while outside_box
                x = bounding_rect_size .* rand(randgen,2) .+ origin(bounding_rect)
                particles[n] = Particle(particle_medium, congruent(particle_shape, x))
                outside_box = !(particles[n] ⊆ box_shape)
            end

            overlapping = false
            for i = 1:(n-1) #compare with previous particles
                if norm(origin(particles[n]) - origin(particles[i])) < separation_ratio * (outer_radius(particles[n]) + outer_radius(particles[i]))
                    overlapping = true
                    break
                end
            end

            num_attempts += 1
            if num_attempts > MAX_ATTEMPTS_TO_FIT_PARTICLE
                error("Tried to place a scatterer $MAX_ATTEMPTS_TO_FIT_PARTICLE times unsuccessfully, just not enough room! You could try increaseing MAX_ATTEMPTS_TO_FIT_PARTICLE")
            end

        end
    end

    return particles
end

function random_particles(particle_medium::P, particle_shape::S,
    box_shape::Shape{T,Dim}, volfrac::Number; seed=Base.Random.make_seed()) where {T,Dim,P<:PhysicalProperties{T,Dim},S<:Shape{T,Dim}}

    N = Int(round(volume(box_shape) / volume(particle_shape) * volfrac))
    return random_particles(particle_medium, particle_shape, box_shape, N; seed=seed)
end
