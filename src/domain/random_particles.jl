MAX_ATTEMPTS_TO_FIT_PARTICLE = 2000

"Generate a non-overlapping Vector of N random particles of radius r that fit inside the shape passed in"
function random_particles(volfrac::Number,a::Number,shape::Shape,seed::Vector{UInt32})
    if volfrac > 0.7854 
        error("Specified volume fraction is larger than optimal packing of circles.")
    end
    N = Int(round(volume(shape) / (π * a^2) * volfrac))
    particles = array_of_particles(N, a)
    random_particles!(particles, N, shape, seed)
    return particles
end

"Generate a non-overlapping Vector of N random particles of radius r that fit inside the shape passed in"
function random_particles(N::Int, a::Number, shape::Shape, seed::Vector{UInt32})
    particles = Particles{typeof(a)}(N, a)
    random_particles!(particles, N, shape, seed)
    return particles
end

"Place N non-overlapping random particles that fit inside the shape passed in. Modifies the particles Vector."
function random_particles!{T}(particles::Vector{Particle{T}}, N::Int, shape::Shape, seed::Vector{UInt32})
    randgen = MersenneTwister(seed)
    # separation distance, the minimum distance between the centres of particles relative to the two radiuses
    separation_ratio = 1.05
    # The factor 2.1radius is seen in several papers (therefore ratio of 1.05)
    println("Generating $N particles with mean radius $(mean_radius(particles)) and total volume of $(volume(particles)) in a shape of volume $(volume(shape)).")
    
    shape_bounding_box = bounding_box(shape)
    topright = shape_bounding_box.topright
    bottomleft = shape_bounding_box.bottomleft
    println("Using a bounding box with volume $(volume(shape_bounding_box))")
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
                error("Tried to place a scatterer $MAX_ATTEMPTS_TO_FIT_PARTICLE times unsuccessfully, just not enough room!")
            end
        end
    end
end