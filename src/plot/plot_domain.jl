@recipe function plot(p::CapsuleParticle)
    @series shape(p.outer)
    @series shape(p.inner)
end

# Plot a vector of particles
@recipe function plot(particles::AbstractParticles)
    for particle in particles
        @series begin
            particle
        end
    end
end

# Plot the shape of a particle
@recipe function plot(particle::AbstractParticle)
    shape(particle)
end


@recipe function plot(shape::Shape{T,2}) where T

    grid --> false
    xguide --> "x"
    yguide --> "y"
    aspect_ratio := 1.0
    label --> ""
    linecolor --> :grey

    x, y = boundary_functions(shape)

    (x, y, 0.0, 1.0)

end

@recipe function plot(shape::Shape{T,3}) where T

    println("Only the (x,z) coordinates of this shape will be drawn")

    grid --> false
    xguide --> "x"
    yguide --> "y"
    aspect_ratio := 1.0
    label --> ""
    linecolor --> :grey

    x, y, z = boundary_functions(shape)

    (x, z, 0.0, 1.0)

end
