const MAX_ATTEMPTS_TO_FIT_PARTICLE = 3000


random_particles(particle_medium::PhysicalProperties{T,Dim}, particle_shape::Shape{T,Dim}; kws...) where {T<:AbstractFloat,Dim} = random_particles(particle_medium, [particle_shape]; kws...)


function random_particles(particle_medium::PhysicalProperties{T,Dim}, particle_shapes::Vector{S};
        num_particles::Int = length(particle_shapes), volume_fraction::T = zero(T),
        box_shape::Shape{T,Dim} = Rectangle(zeros(T,2), T(10)*sum(outer_radius.(particle_shapes)), T(10)*sum(outer_radius.(particle_shapes))),
        current_particles::Vector{AbstractParticle{T,Dim}} = AbstractParticle{T,Dim}[],
        kws...
) where {T<:AbstractFloat,Dim, S<:Shape{T,Dim}}

    particles = if volume_fraction == zero(T)
        num_each = Int(round(num_particles / length(particle_shapes)))
        for s in particle_shapes
            add_particles = random_particles(particle_medium, s, box_shape, num_each; current_particles = current_particles, kws...)
            current_particles = [current_particles; add_particles]
        end
        current_particles
    else
        vol_each = volume_fraction / length(particle_shapes)
        for s in particle_shapes
            add_particles = random_particles(particle_medium, s, box_shape, vol_each; current_particles = current_particles, kws...)
            current_particles = [current_particles; add_particles]
        end
        current_particles
    end

    return particles
end

# random_particles(particle_medium, particle_shapes)

#
# a = 1
# function f1(; a = [])
#     if volume_fraction != zero(T)
#         volume_fraction_each = volume_fraction / length(particle_shapes)
#         for s in particle_shapes
#             b = a
#             a = [a; b]
#             println(a)
#         end
#         println(a)
#     end
#     println(a)
# end
# f1(; a = [2])
# a

# function random_particles(particle_medium::PhysicalProperties{T,Dim}, particle_shape::Shape{T,Dim};
#         box_shape::Shape{T,Dim} = Rectangle(zeros(T,2), T(10)*outer_radius(particle_shape), T(10)*outer_radius(particle_shape)),
#         num_particles::Int = 0, volume_fraction::T = zero(T), kws...) where {T<:AbstractFloat,Dim}
#
#     if volume_fraction == zero(T)
#         if num_particles == 0 num_particles = 5 end
#         random_particles(particle_medium, particle_shape,
#             box_shape, num_particles; kws...)
#     else
#         random_particles(particle_medium, particle_shape,
#             box_shape, volume_fraction; kws...)
#     end
#
# end

function random_particles(particle_medium::P, particle_shape::S,
    box_shape::Shape{T,Dim}, volfrac::AbstractFloat; kws...
) where {T,Dim,P<:PhysicalProperties{T,Dim},S<:Shape{T,Dim}}

    N = Int(round(volume(box_shape) / volume(particle_shape) * volfrac))
    return random_particles(particle_medium, particle_shape, box_shape, N; kws...)
end

"""
    random_particles(particle_medium, particle_shapes::Vector{Shape}, box_shape, volume_fraction::Number;
        seed=Random.make_seed())
    random_particles(particle_medium, particle_shape::Shape, box_shape, volume_fraction::Number;
        seed=Random.make_seed())
    random_particles(particle_medium, particle_shape::Shape, box_shape, N::Integer;
        seed=Random.make_seed())

Generate `N` random particles that fit inside `box_shape` (or fill with `volume_fraction`)

Specify seed to make output deterministic. Algorithm places particles unifomly randomly inside
the bounding rectangle of `box_shape` and discards particle if it overlaps (based on outer radius)
or does not lies completely in box box.

When passing particle_shapes::Vector{Shape} we assume each element is equally likely to occur. Repeating the same shape will lead to it being placed more often.
"""
function random_particles(particle_medium::P, particle_shape::S, box_shape::Shape{T,Dim}, N::Integer;
        seed=Random.make_seed(),
        separation_ratio::T = T(1.005), # Min distance between particle centres relative to their outer radiuses.
        current_particles::Vector{AbstractParticle{T,Dim}} = AbstractParticle{T,Dim}[] # Particles already present.
) where {T<:AbstractFloat,Dim,P<:PhysicalProperties{T,Dim},S<:Shape{T,Dim}}

    # Check volume fraction is not impossible
    volfrac = N * volume(particle_shape) / volume(box_shape)
    max_packing = if length(current_particles) > 0
        0.7854 - sum(volume.(current_particles))/volume(box_shape)
    else 0.7854
    end
    if volfrac  > max_packing
        error("Specified volume fraction is larger than optimal packing of circles.")
    end

    # Create pseudorandom device with specified seed
    randgen = MersenneTwister(seed)

    bounding_rect = bounding_rectangle(box_shape)
    bounding_rect_size = SVector(bounding_rect.width, bounding_rect.height)

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
    L = length(current_particles)
    particles = Vector{Particle{T,Dim,P,S}}(undef, N + L)
    particles[1:L] = current_particles

    for n = L+1:L+N

        num_attempts = 0
        overlapping = true
        while overlapping

            outside_box = true
            while outside_box
                x = bounding_rect_size .* (1 .- 2 .* rand(randgen,T,2)) .+ origin(bounding_rect)
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
                error("Tried to place a scatterer $MAX_ATTEMPTS_TO_FIT_PARTICLE times unsuccessfully. You could try increasing MAX_ATTEMPTS_TO_FIT_PARTICLE")
            end

        end
    end

    return particles[L+1:L+N]
end
