@recipe function plot(particle::Particle) @series particle.shape end

#Currently plot(Vector{Particle}) does not work
@recipe function plot(particles::Vector{P}) where P <: Particle
    # for p in particles @series identity(p) end
    for i=1:length(particles) @series particles[i] end

    # @series particles[1]
end

@recipe function plot(shape::Shape)
    grid --> false
    legend --> nothing
    xlab --> "x"
    ylab --> "y"
    aspect_ratio := 1.0
    label --> name(shape)
    fill --> (0, :transparent)
    linecolor --> :grey

    x, y = boundary_functions(shape)

    (x, y, 0, 1)
end
