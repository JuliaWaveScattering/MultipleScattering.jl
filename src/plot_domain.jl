@recipe function plot(particles::Vector{Particle})
    # for p in particles @series identity(p) end
    # for i=1:length(particles) @series particles[i] end

    @series particles[1]
end

@recipe function plot(particle::Particle)
    grid --> false
    legend --> nothing
    xlab --> "x"
    ylab --> "y"
    aspect_ratio := 1.0
    fillalpha = 0.0
    # fillalpha = 0.7/(1.0 + abs(particle.ρ*particle.c)) # darker fill the larger the impendence
    fill --> (0, fillalpha, :grey)
    linecolor --> :grey
    linealpha --> 0.7

    x = t -> particle.r*cos(t) + particle.x[1]
    y = t -> particle.r*sin(t) + particle.x[2]

    (x, y, 0, 2π)
end

@recipe function plot(shape::Shape)
    grid --> false
    legend --> nothing
    xlab --> "x"
    ylab --> "y"
    aspect_ratio := 1.0
    label --> name(shape)
    fill --> (0, :transparent)
    line --> (1, :red)

    x, y = boundary_functions(shape)

    (x, y, 0, 1)
end
