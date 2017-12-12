@recipe function plot(particles::Vector{Particle})
    grid --> false
    legend --> nothing
    xlab --> "x"
    ylab --> "y"
    aspect_ratio := 1.0
    fill --> (0, :grey)
    line --> 0

    x = map(p -> (t -> p.r*cos(t) + p.x[1]), particles)
    y = map(p -> (t -> p.r*sin(t) + p.x[2]), particles)

    (x, y, 0, 2π)
end

@recipe function plot(particle::Particle)
    grid --> false
    legend --> nothing
    xlab --> "x"
    ylab --> "y"
    aspect_ratio := 1.0
    fillalpha = 0.7/(1.0 + abs(particle.ρ*particle.c)) # darker fill the larger the impendence
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
