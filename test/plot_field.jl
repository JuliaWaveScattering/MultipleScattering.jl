"Plot just the particles and source"
function plot_field(simres::SimulationResult; seriestype = :line, field_apply = real)
    if seriestype == :line
        labs = [ "real x=$x" for x in simres.x];
        labs = [labs ; [ "imag x=$x" for x in simres.x]]
        ωt = getfield(simres, fieldnames(simres)[end])
        plot(transpose(ωt), transpose([real.(field(simres)); imag.(field(simres))])
            , labels = reshape(labs,1,length(labs))
        )

    elseif seriestype == :contour || seriestype == :field
        ind_ωt = 1 # choose which time or frequency

        x_pixels = union([x[1] for x in simres.x])
        y_pixels = union([x[2] for x in simres.x])

        # Make the vector field(field_sim)[:,ind_ωt] into a matrix for heatmap
        field_mat = transpose(
            reshape(field(simres)[:,ind_ωt], (length(x_pixels), length(y_pixels)))
        )
        # linetype --> :contour
        # fill --> true
        # aspect_ratio := 1.0
        # fillcolor --> :pu_or
        # title --> "Field at ω=$ω"

        plot(x_pixels, y_pixels, field_apply.(field_mat)
            , fill=true, aspect_ratio = 1.0, fillcolor = :pu_or, seriestype = :contour)
    else error("Unknown linestyle = $linestyle for $res")
    end
end
