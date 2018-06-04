include("plot_domain.jl")
# include("plot_moments.jl")

"Plot just the particles and source"
@recipe function plot(simres::SimulationResult; linetype = :line, field_apply = real)

    if linetype == :line
        @series begin
            labs = [ "real x=$x" for x in simres.x];
            labs = [labs ; [ "imag x=$x" for x in simres.x]]
            ωt = getfield(simres, fieldnames(simres)[end])

            labels --> reshape(labs,1,length(labs))
            (transpose(ωt), transpose([real.(field(simres)); imag.(field(simres))]))
        end

    elseif linetype == :contour || linetype == :field
        @series begin
            ind_ωt = 1 # choose which time or frequency

            x_pixels = union([x[1] for x in simres.x])
            y_pixels = union([x[2] for x in simres.x])

            # Make the vector field(field_sim)[:,ind_ωt] into a matrix for heatmap
            field_mat = transpose(
                reshape(field(simres)[:,ind_ωt], (length(x_pixels), length(y_pixels)))
            )
            linetype --> :contour
            fill --> true
            aspect_ratio := 1.0
            fillcolor --> :pu_or
            # title --> "Field at ω=$ω"

            (x_pixels, y_pixels, field_apply.(field_mat))
        end
    else error("Unknown linestyle = $linestyle for $res")
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
        linetype --> :contour
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


# "Plot the response across all wavenumbers"
# @recipe function plot(simulation::FrequencySimulation)
#     label --> ["real" "imaginary"]
#     xlabel --> "Wavenumber (k)"
#     ylabel --> "Response"
#     grid --> false
#     title --> "Response from particles of radius $(signif(simulation.particles[1].r,2)) contained in a $(lowercase(name(simulation.shape)))\n with volfrac=$(signif(calculate_volfrac(simulation),2)) measured at ($(simulation.listener_positions[1,1]), $(simulation.listener_positions[2,1]))"
#
#     (simulation.k_arr, [real(simulation.response) imag(simulation.response)])
# end
