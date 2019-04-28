# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(simres::FrequencySimulationResult, ω::AbstractFloat;
        time = 0.0, # does not except "t" as name of variable..
        x_indices = axes(simres.x,1),
        ω_index = findmin(abs.(getfield(simres, 3) .- ω))[2],
        field_apply = real, seriestype = :surface)

    x = [x[1] for x in simres.x[x_indices]]
    y = [x[2] for x in simres.x[x_indices]]
    ω = getfield(simres, 3)[ω_index]

    phase = exp(-im*ω*time)

    color --> :pu_or
    title --> "Field for ω = $ω"
    seriestype --> seriestype
    aspect_ratio --> 1.0

    if seriestype == :contour
        # We should really check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true
        x, y, field_apply.(phase.*transpose(reshape(field(simres)[x_indices,ω_index],n_y,n_x)))

    else
        (x, y, field_apply.(phase.*field(simres)[x_indices,ω_index]))
    end

end

# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(timres::TimeSimulationResult, t::AbstractFloat;
        x_indices = axes(timres.x,1),
        t_index = findmin(abs.(getfield(timres, 3) .- t))[2],
        field_apply = real, seriestype = :surface)

    x = [x[1] for x in timres.x[x_indices]]
    y = [x[2] for x in timres.x[x_indices]]
    t = getfield(timres, 3)[t_index]

    color --> :pu_or
    title --> "Field for time = $t"
    seriestype --> seriestype
    aspect_ratio --> 1.0

    if seriestype == :contour

        # We should really check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true
        x, y, field_apply.(transpose(reshape(field(timres)[x_indices,t_index],n_y,n_x)))

    else

        (x, y, field_apply.(field(timres)[x_indices,t_index]))

    end

end
