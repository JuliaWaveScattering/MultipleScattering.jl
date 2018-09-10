include("plot_domain.jl")

# Plot the result across angular frequency for a specific position (x)
@recipe function plot(simres::SimulationResult;
        x = simres.x,
        x_indices = [findmin(norm(z - y) for z in simres.x)[2] for y in x],
        ω_indices = Colon(), apply = real)

    for x_ind in x_indices

        apply_field = apply.(field(simres)[x_ind, ω_indices])
        xlab = ((typeof(simres) <: FrequencySimulationResult) ? "ω" : "t")

        @series begin
            label --> "$apply x=$(simres.x[x_ind])"
            xlabel --> xlab
            (transpose(getfield(simres, 3)[ω_indices]), apply_field)
        end
    end
end

# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(simres::SimulationResult, ω_or_t::AbstractFloat;
        x_indices = indices(simres.x,1),
        ω_or_t_index = findmin(abs.(getfield(simres, 3) .- ω_or_t))[2],
        field_apply = real, seriestype = :surface)

    x = [x[1] for x in simres.x[x_indices]]
    y = [x[2] for x in simres.x[x_indices]]
    ω_or_t = getfield(simres, 3)[ω_or_t_index]

    fillcolor --> :pu_or
    title --> "Field at $(fieldnames(simres)[end])=$ω_or_t"
    seriestype --> seriestype
    aspect_ratio --> 1.0

    if seriestype == :contour

        # We should really check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true
        x, y, field_apply.(transpose(reshape(field(simres)[x_indices,ω_or_t_index],n_y,n_x)))

    else

        (x, y, field_apply.(field(simres)[x_indices,ω_or_t_index]))

    end

end

"Plot just the particles"
@recipe function plot(sim::FrequencySimulation; bounds = :none)

    println("Plotting a simulation on its own")

    @series begin
        if bounds != :none
            # bounds = bounding_rectangle(sim.particles)
            xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
            ylims --> (bottomleft(bounds)[2], topright(bounds)[2])
        end

        sim.particles
    end

end

"Plot the field for a particular wavenumber"
@recipe function plot(sim::FrequencySimulation{T}, ω::T; res=10, xres=res, yres=res,
                         field_apply=real, bounds = :auto, drawparticles=false) where {T}

    # If user wants us to, generate bounding rectangle around particles
    if bounds == :auto
        if isempty(sim.particles)
            warn("What region to plot? Use keyword bounds = Rectangle")
            bounds = Rectangle([-one(T),-one(T)],[one(T),one(T)])
        else
            bounds = bounding_rectangle(sim.particles)
        end
    end

    # If user has not set xlims and ylims, set them to the rectangle
    xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
    ylims --> (bottomleft(bounds)[2], topright(bounds)[2])

    # Incase the user did set the xlims and ylims, generate a new bounding
    # rectangle with them
    p_xlims = plotattributes[:xlims]
    p_ylims = plotattributes[:ylims]
    bounds = Rectangle((T(p_xlims[1]),T(p_ylims[1])), (T(p_xlims[2]),T(p_ylims[2])))

    @series begin

        field_sim = run(sim, bounds, [ω]; xres=xres, yres=yres)

        xy_mat = reshape(field_sim.x, (xres+1, yres+1))

        x_pixels = [x[1] for x in xy_mat[:,1]]
        y_pixels = [x[2] for x in xy_mat[1,:]]

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field(field_sim), (xres+1, yres+1)))
        seriestype --> :contour
        fill --> true
        grid --> false
        aspect_ratio := 1.0
        color --> :pu_or
        title --> "Field at ω=$ω"

        (x_pixels, y_pixels, field_apply.(response_mat))
    end

    if drawparticles
        @series begin
            sim.particles
        end
    end

end
