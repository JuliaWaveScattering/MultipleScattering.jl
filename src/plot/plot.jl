include("plot_domain.jl")
include("plot_field.jl")
include("plot_moments.jl")

# Plot the result across angular frequency for a specific position (x)
@recipe function plot(simres::SimulationResult;
        x = simres.x,
        x_indices = [findmin([norm(z - y) for z in simres.x])[2] for y in x],
        ω_indices = Colon(), field_apply = real)

    for x_ind in x_indices

        fs = field_apply.(field(simres)[x_ind, ω_indices])
        xguide = ((typeof(simres) <: FrequencySimulationResult) ? "ω" : "t")

        @series begin
            label --> "$field_apply x=$(simres.x[x_ind])"
            xguide --> xguide
            (getfield(simres, 3)[ω_indices], fs)
        end
    end
end

"Plot just the particles"
@recipe function plot(sim::FrequencySimulation; bounds = :none)

    # println("Plotting a simulation on its own")

    @series begin
        if bounds != :none
            # bounds = bounding_box(sim.particles)
            xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
            ylims --> (bottomleft(bounds)[2], topright(bounds)[2])
        end

        sim.particles
    end

end


"Plot the field for a particular wavenumber"
@recipe function plot(sim::FrequencySimulation{T,3}, ω::T;
        res=10, xres=res, yres=res,
        y = :auto,
        field_apply=real,
        region_shape = :auto,
        bounds = :auto,
        exclude_region = EmptyShape{T,3}(),
        drawparticles=false) where {T}

    # If user wants us to, generate bounding rectangle around particles
    region_shape = (region_shape != :auto) ? region_shape :
        if isempty(sim.particles)
            if bounds == :auto
                @warn "What region to plot? For example, use keyword bounds = Box([[-1.0,-1.0],[1.0,1.0]])"
                Box([[-one(T),-one(T)],[one(T),one(T)]])
            else bounds
            end
        else
            region_shape = bounding_box(sim.particles)
        end

    bounds = bounding_box(region_shape)

    # If user has not set xlims and ylims, set them to the rectangle
    corners = box_corners(bounds)
    xlims --> (corners[1][1], corners[end][1])
    ylims --> (corners[1][end], corners[end][end])

    if y == :auto
        y = corners[1][2]
    end

    # Incase the user did set the xlims and ylims, generate a new bounding box with them
    # p_xlims = plotattributes[:xlims]
    # p_ylims = plotattributes[:ylims]
    # bounds = Box([[T(p_xlims[1]),T(p_ylims[1])], [T(p_xlims[2]),T(p_ylims[2])]])
    #
    # region_shape = (bounds ⊆ region_shape) ? bounds : region_shape

    field_sim = run(sim, region_shape, [ω]; y=y, xres=xres, zres=yres, exclude_region=exclude_region)
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
@recipe function plot(sim::FrequencySimulation{T,2}, ω::T;
        res=10, xres=res, yres=res,
        field_apply=real,
        region_shape = :auto,
        bounds = :auto,
        exclude_region = EmptyShape{T,2}(),
        drawparticles=false) where {T}

    # If user wants us to, generate bounding rectangle around particles
    region_shape = (region_shape != :auto) ? region_shape :
        if isempty(sim.particles)
            if bounds == :auto
                @warn "What region to plot? For example, use keyword bounds = Box([[-1.0,-1.0],[1.0,1.0]])"
                Box([[-one(T),-one(T)],[one(T),one(T)]])
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
    bounds = Box([[T(p_xlims[1]),T(p_ylims[1])], [T(p_xlims[2]),T(p_ylims[2])]])

    region_shape = (bounds ⊆ region_shape) ? bounds : region_shape

    field_sim = run(sim, region_shape, [ω]; xres=xres, yres=yres, zres=zres, exclude_region=exclude_region)
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
@recipe function plot(s::Source, ω)
    (FrequencySimulation(s), ω)
end
