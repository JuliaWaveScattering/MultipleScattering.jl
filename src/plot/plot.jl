include("plot_domain.jl")
# include("plot_moments.jl")

# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(simres::FrequencySimulationResult, x_indices::Union{UnitRange{Int},Colon,AbstractVector{Int}}, ω_index::Int; field_apply = real, seriestype=:surface)

    x = [x[1] for x in simres.x[x_indices]]
    y = [x[2] for x in simres.x[x_indices]]

    if seriestype == :contour

        # We should really check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true
        fillcolor --> :pu_or

        x, y, field_apply.(reshape(field(simres)[x_indices,ω_index],n_y,n_x))

    else

        seriestype := :surface
        fillcolor --> :pu_or
        title --> "Field at ω=$(simres.ω[ω_index])"

        (x, y, field_apply.(field(simres)[x_indices,ω_index]))

    end

end

# Plot the result across angular frequency for a specific position (x)
@recipe function plot(simres::FrequencySimulationResult, x_index::Int, ω_indices::Union{UnitRange{Int},Colon,AbstractVector{Int}}=Colon(); field_apply = real)

    x = simres.x[x_index]
    complex_field = field(simres)[x_index, ω_indices]

    @series begin
        label --> "Real x=$x"
        (transpose(simres.ω[ω_indices]), real.(complex_field))
    end

    @series begin
        label --> "Imag x=$x"
        (transpose(simres.ω[ω_indices]), imag.(complex_field))
    end

end

"Plot just the particles"
@recipe function plot(sim::FrequencySimulation; bounds = :auto)

    println("Plotting a simulation on its own")

    if bounds == :auto
        bounds = bounding_rectangle(sim.particles)
    end

    @series begin
        xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
        ylims --> (bottomleft(bounds)[2], topright(bounds)[2])
        sim.particles
    end

end

"Plot the field for a particular wavenumber"
@recipe function plot(sim::FrequencySimulation, ω::Number; res=10, xres=res, yres=res,
                         field_apply=real, bounds = :auto, drawparticles=false)

    # If user wants us to, generate bounding rectangle around particles
    if bounds == :auto
        bounding_rect = bounding_rectangle(sim.particles)
    end

    # If user has not set xlims and ylims, set them to the rectangle
    xlims --> (bottomleft(bounding_rect)[1], topright(bounding_rect)[1])
    ylims --> (bottomleft(bounding_rect)[2], topright(bounding_rect)[2])

    # Incase the user did set the xlims and ylims, generate a new bounding
    # rectangle with them
    p_xlims = plotattributes[:xlims]
    p_ylims = plotattributes[:ylims]
    bounding_rect = Rectangle([p_xlims[1],p_ylims[1]], [p_xlims[2],p_ylims[2]])

    @series begin

        field_sim = run(sim, bounding_rect, [ω]; xres=xres, yres=yres)

        xy_mat = reshape(field_sim.x, (xres+1, yres+1))

        x_pixels = [x[1] for x in xy_mat[:,1]]
        y_pixels = [x[2] for x in xy_mat[1,:]]

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field(field_sim), (xres+1, yres+1)))
        seriestype --> :contour
        fill --> true
        aspect_ratio := 1.0
        fillcolor --> :pu_or
        title --> "Field at ω=$ω"

        (x_pixels, y_pixels, field_apply.(response_mat))
    end
    if drawparticles
        @series begin
            sim.particles
        end
    end
end
