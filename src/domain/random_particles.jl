const MAX_ATTEMPTS_TO_FIT_PARTICLE = 2000

"Generate a non-overlapping Vector of N random particles of radius r that fit inside the shape passed in"
function random_particles{T}(volfrac::Number, a::T, shape::Shape;
        c=one(Complex{T}), ρ=zero(T),
        seed::Vector{UInt32}=Base.Random.make_seed()
    )
    if volfrac > 0.7854
        error("Specified volume fraction is larger than optimal packing of circles.")
    end
    N = Int(round(volume(shape) / (π * a^2) * volfrac))
    particles = [Particle{typeof(a)}(zeros(T, 2), a, c, ρ) for i=1:N]
    random_particles!(particles, shape; seed=seed)
    return particles
end

"Generate a non-overlapping Vector of N random particles of radius r that fit inside the shape passed in"
function random_particles{T}(N::Int, a::T, shape::Shape = Rectangle(0.1, a, N);
        c=one(Complex{T}), ρ=zero(T),
        seed::Vector{UInt32}=Base.Random.make_seed()
    )
    # particles = Particles{typeof(a)}(N, a)
    particles = [Particle{typeof(a)}(zeros(T, 2), a, c, ρ) for i=1:N]
    random_particles!(particles, shape; seed=seed)
    return particles
end

"Place N non-overlapping random particles that fit inside the shape passed in. Modifies the particles Vector."
function random_particles!{T}(particles::Vector{Particle{T}}, shape::Shape;
        N::Int=length(particles),
        seed::Vector{UInt32}=Base.Random.make_seed()
    )
    randgen = MersenneTwister(seed)
    # separation distance, the minimum distance between the centres of particles relative to the two radiuses
    separation_ratio = 1.05
    # The factor 2.1radius is seen in several papers (therefore ratio of 1.05)

    shape_bounding_box = bounding_box(shape)
    topright = shape_bounding_box.topright
    bottomleft = shape_bounding_box.bottomleft
    
    @printf("Generating %d randomly positioned particles with mean radius of %0.5g. Total volume of particles is %0.5g in a %s of volume %0.5g (volume fraction %0.5g). Using a bounding box with volume %0.5g. \n", N, mean_radius(particles), volume(particles), name(shape), volume(shape), volume(particles)/volume(shape), volume(shape_bounding_box))
    
    for n = 1:N

        num_attempts = 0
        overlapping = true
        while overlapping
            # Generate a random position inside the bounding box from the seeded random number generator, then test if it is the shape, if not try again
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
