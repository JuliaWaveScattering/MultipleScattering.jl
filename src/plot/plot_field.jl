
"Plot the field for a particular wavenumber"
@recipe function plot(sim::FrequencySimulation{3}, ω::Number;
        resolution = 10, res = resolution, xres=res, yres=res,
        y = :auto,
        field_apply=real,
        region_shape = :auto,
        bounds = :auto,
        exclude_region = EmptyShape{3}(),
        drawparticles=false)

    # If user wants us to, generate bounding rectangle around particles
    region_shape = (region_shape != :auto) ? region_shape :
        if isempty(sim.particles)
            if bounds == :auto
                @warn "What region to plot? For example, use keyword bounds = Box([[-1.0,-1.0],[1.0,1.0]])"
                Box([[-1,-1],[1,1]])
            else bounds
            end
        else
            region_shape = bounding_box(sim.particles)
        end

    bounds = bounding_box(region_shape)

    # If user has not set xlims and ylims, set them to the rectangle
    cs = corners(bounds)
    xlims --> (cs[1][1], cs[end][1])
    ylims --> (cs[1][end], cs[end][end])

    if y == :auto
        y = cs[1][2]
    end

    # Incase the user did set the xlims and ylims, generate a new bounding box with them
    # p_xlims = plotattributes[:xlims]
    # p_ylims = plotattributes[:ylims]
    # bounds = Box([[T(p_xlims[1]),T(p_ylims[1])], [T(p_xlims[2]),T(p_ylims[2])]])
    #
    # region_shape = (bounds ⊆ region_shape) ? bounds : region_shape

    field_sim = run(sim, region_shape, [ω];
        y=y, xres=xres, zres=yres,
        exclude_region=exclude_region
    )

    xy_mat = reshape(field_sim.x, (xres+1, yres+1))
    x_pixels = [x[1] for x in xy_mat[:,1]]
    y_pixels = [x[end] for x in xy_mat[1,:]]

    @series begin

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field(field_sim), (xres+1, yres+1)))
        seriestype --> :contour
        fill --> true
        grid --> false
        aspect_ratio := 1.0
        seriescolor --> :balance
        title --> "Field at ω=$ω"

        (x_pixels, y_pixels, field_apply.(response_mat))
    end

    if drawparticles
        @series begin
            sim.particles
        end
    end

end


"Plot the field for a particular wavenumber"
@recipe function plot(sim::FrequencySimulation{2}, ω::Number;
        resolution = 10, res = resolution, xres=res, yres=res,
        field_apply=real,
        region_shape = :auto,
        bounds = :auto,
        exclude_region = EmptyShape{2}(),
        drawparticles=true)

    # If user wants us to, generate bounding rectangle around particles
    region_shape = (region_shape != :auto) ? region_shape :
        if isempty(sim.particles)
            if bounds == :auto
                @warn "What region to plot? For example, use keyword bounds = Box([[-1.0,-1.0],[1.0,1.0]])"
                Box([[-1.0,-1.0],[1.0,1.0]])
            else bounds
            end
        else
            region_shape = bounding_box(sim.particles)
        end

    bounds = bounding_box(region_shape)
    # If user has not set xlims and ylims, set them to the rectangle
    xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
    ylims --> (bottomleft(bounds)[2], topright(bounds)[2])

    # Incase the user did set the xlims and ylims, generate a new bounding
    # rectangle with them
    p_xlims = plotattributes[:xlims]
    p_ylims = plotattributes[:ylims]
    bounds = Box([[p_xlims[1],p_ylims[1]], [p_xlims[2],p_ylims[2]]])

    region_shape = (bounds ⊆ region_shape) ? bounds : region_shape

    field_sim = run(sim, region_shape, [ω]; xres=xres, yres=yres, zres=yres, exclude_region=exclude_region)
    xy_mat = reshape(field_sim.x, (xres+1, yres+1))
    x_pixels = [x[1] for x in xy_mat[:,1]]
    y_pixels = [x[2] for x in xy_mat[1,:]]

    @series begin

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field(field_sim), (xres+1, yres+1)))
        seriestype --> :contour
        fill --> true
        grid --> false
        aspect_ratio := 1.0
        seriescolor --> :balance
        title --> "Field at ω=$ω"

        (x_pixels, y_pixels, field_apply.(response_mat))
    end

    if drawparticles
        @series begin
            sim.particles
        end
    end

end

"Plot the source field for a particular wavenumber"
@recipe function plot(s::RegularSource, ω)
    (FrequencySimulation(s), ω)
end

# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(simres::FrequencySimulationResult, ω::AbstractFloat;
        phase_time = 0.0, # does not except "t" as name of variable..
        x_indices = axes(simres.x,1),
        ω_index = findmin(abs.(getfield(simres, 3) .- ω))[2],
        field_apply = real, seriestype = :surface,
        region_shape = :empty)

    x = [x[1] for x in simres.x[x_indices]]

    # y will actually be z for 3D...
    y = [x[end] for x in simres.x[x_indices]]
    ω = getfield(simres, 3)[ω_index]

    phase = exp(-im*ω*phase_time)

    seriestype --> seriestype
    seriescolor --> :balance
    title --> "Field for ω = $ω"
    aspect_ratio --> 1.0

    if seriestype != :surface
        # We could check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true

        if region_shape != :empty
            bounds = bounding_box(region_shape)

            # If user has not set xlims and ylims, set them to the rectangle
            xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
            ylims --> (bottomleft(bounds)[2], topright(bounds)[2])
        else
            xlims --> (minimum(x), maximum(x))
            ylims --> (minimum(y), maximum(y))
        end

        x, y, field_apply.(phase.*transpose(reshape(field(simres)[x_indices,ω_index],n_x,n_y)))

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

    seriescolor --> :balance
    title --> "Field for time = $t"
    seriestype --> seriestype
    aspect_ratio --> 1.0

    if seriestype != :surface

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
