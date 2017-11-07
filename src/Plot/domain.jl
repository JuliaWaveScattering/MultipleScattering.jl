
# "Draw each particle as a circle, ensuring it is true to size and shape"
function plot_particles{T}(particles::Vector{Particle{T}}; newplot=true)
    # Make vector of function pointers for each particle, which will draw a circle centred at said particle
    P = length(particles)
    if newplot plot() end

    if P==0
        warn("There are no particles to draw")
        plot!()
    else
        x_particles_boundaries = Array{Function, 1}(P)
        y_particles_boundaries = Array{Function, 1}(P)
        for p=1:P
            x_particles_boundaries[p] = t->particles[p].r * sin(t) + particles[p].x[1]
            y_particles_boundaries[p] = t->particles[p].r * cos(t) + particles[p].x[2]
        end
        # Add all the particles to the plot as paramteric plots
        plot!(x_particles_boundaries, y_particles_boundaries, 0, 2π, line=0, fill=(0, :grey), aspect_ratio=1.0)
    end
    plot!(grid=false, legend=nothing, xlab="x", ylab="y")
end

# "Plots the particles from a specific model. Wrapper using model instead of vector of particles to plot the particles"
function plot_particles{T}(model::FrequencyModel{T}; newplot=true)
    plot_particles(model.particles;newplot=newplot)
end

# "Draw the listener as a pacman with a radius the mean of all the particles"
function plot_listeners{T}(model::FrequencyModel{T}; newplot=true)
    if newplot plot() end

    if length(model.particles)==0
        warn("There are no particles, the listener will be drawn with a radius of 0.5")
        r = .5
    else
        r = mean_radius(model.particles)/2.0
    end

    for i=1:size(model.listener_positions, 2)
        x_listener_boundary(t) = r * cos(t) + model.listener_positions[1, i]
        y_listener_boundary(t) = r * sin(t) + model.listener_positions[2, i]
        plot!(x_listener_boundary, y_listener_boundary, -π/3.,π/3., line=0, fill=(0, :blue), aspect_ratio=1.0)
    end
    plot!()
end

function plot_shape(shape::Shape)
    warn("plot_shape not implemented")
end

# "Plot all the elements of the domain"
function plot_domain{T}(model::FrequencyModel{T}; newplot=true, particles_flag=true, shape_flag=true, listeners_flag=true)
    if newplot
        plot()
    end

    if particles_flag
        plot_particles(model; newplot=false)
    end

    if shape_flag
        plot_shape(model.shape)
    end

    if listeners_flag
        plot_listeners(model; newplot=false)
    end

    plot!()
end
