@recipe function plot(particle::AbstractParticle) @series shape(particle) end

@recipe function plot(p::CapsuleParticle)
    @series shape(p.outer)
    @series shape(p.inner)
end

#Currently plot(Vector{Particle}) does not work
@recipe function plot(particles::Vector{P}) where P <: AbstractParticle
    # for p in particles @series identity(p) end
    for i=1:length(particles) @series particles[i] end
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
