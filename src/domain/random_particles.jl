const MAX_ATTEMPTS_TO_FIT_PARTICLE = 2000

"""
Generate a non-overlapping Vector of N random particles of radius r that fit 
inside the shape passed in
"""
function random_particles{T}(volfrac::Number, a::T, shape::Shape;
        c=one(Complex{T}), ρ=zero(T),
        seed=Base.Random.make_seed()
    )
    if volfrac > 0.7854
        error("Specified volume fraction is larger than optimal packing of circles.")
    end
    N = Int(round(volume(shape) / (π * a^2) * volfrac))
    return random_particles(N, a, shape;c=c, ρ=ρ, seed=seed)
end

"""
Generate a non-overlapping Vector of N random particles of radius r that fit 
inside the shape passed in
"""
function random_particles{T}(N::Int, a::T, shape::Shape = Rectangle(0.1, a, N);
        c=one(Complex{T}), ρ=zero(T),
        seed=Base.Random.make_seed()
    )
    particles = [Particle{typeof(a)}(zeros(T, 2), a, c, ρ) for i=1:N]
    random_particles!(particles, shape; seed=seed)
    return particles
end

"""
Place N non-overlapping random particles that fit inside the shape passed in. 
Modifies the particles Vector.
"""
function random_particles!{T}(particles::Vector{Particle{T}}, shape::Shape;
        N::Int=length(particles),
        seed=Base.Random.make_seed()
    )
    randgen = MersenneTwister(seed)
    # Min distance between the centre of particles relative to their radiuses
    separation_ratio = 1.05
    # The factor 2.1radius is seen in several papers (therefore ratio of 1.05)

    bounds = bounding_box(shape)
    topright = bounds.topright
    bottomleft = bounds.bottomleft
    
    @printf("""\n
        Generating %d randomly positioned particles
        Mean radius: %0.5g
        Total particle volume: %0.5g
        Inside %s of volume: %0.5g 
        Particle volume fraction: %0.5g 
        Bounding box volume: %0.5g
        """, N, mean_radius(particles), volume(particles), name(shape), 
        volume(shape), volume(particles)/volume(shape), volume(bounds)
    )

    for n = 1:N

        num_attempts = 0
        overlapping = true
        while overlapping
            # Generate a random position inside the bounding box from the seeded
            # random number generator, accept it if in shape, if not try again
            outside = true
            while outside
                particles[n].x .= (topright .- bottomleft) .* rand(randgen,2) .+ bottomleft
                outside = !inside(shape,particles[n])
            end
            # println("Trying to place particle $n at $(particles[n].x)")
            overlapping = false
            for i = 1:(n-1) #compare with previously initialised particles
                # println("Placing particle $n, checking particle $i, separation is $(norm(particles[n].x - particles[i].x))")
                if norm(particles[n].x - particles[i].x) < separation_ratio * (particles[n].r + particles[i].r)
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
end
