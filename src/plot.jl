
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
    fill --> (0, :grey)
    line --> 0

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

"""
Build a 'field model' with lots of listeners using the same domain as model 
you pass in. This 'field model' can then be used to plot the whole field for 
this wavenumber.
"""
function build_field_model{T}(model::FrequencyModel{T},k::T;res=10,xres=res,yres=res)
    # Create deep copy of model so that we can add lots of new listener positions and rerun the model
    field_model = deepcopy(model)

    # Build the listeners or pixels
    bounds = bounding_box(model.shape)
    box_size = bounds.topright - bounds.bottomleft
    box_width = box_size[1]
    box_height = box_size[2]

    # Build up the pixels and all the framework for the plotting
    num_pixels = (xres+1)*(yres+1)
    listener_positions = Matrix{T}(2,num_pixels)

    #Size of the step in x and y direction
    step_size = [box_width / xres, box_height / yres]

    iterator = 1
    for j=0:yres
        for i=0:xres
            listener_positions[:,iterator] = bounds.bottomleft + step_size.*[i,j]
            iterator += 1
        end
    end

    field_model.listener_positions = listener_positions
    field_model.response = Matrix{T}(1, num_pixels)
    generate_responses!(field_model, [k])

    return field_model
end

"Plot the field for a particular wavenumber"
@recipe function plot{T}(model::FrequencyModel{T},k::T;res=10, xres=res, yres=res, resp_fnc=real)

    @series begin
        field_model = build_field_model(model, k; xres=xres, yres=yres)

        # Build the listeners or pixels
        bounds = bounding_box(field_model.shape)

        # For this we sample at the centre of each pixel
        x_pixels = linspace(bounds.bottomleft[1], bounds.topright[1], xres+1)
        y_pixels = linspace(bounds.bottomleft[2], bounds.topright[2], yres+1)

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field_model.response, (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        # fillcolor=fillcolor

        (x_pixels, y_pixels, resp_fnc(response_mat))
    end

    @series begin
        model.shape
    end

    for i=1:length(model.particles) @series model.particles[i] end

    @series begin
        line --> 0
        fill --> (0, :blue)
        legend --> false
        grid --> false
        colorbar --> true
        aspect_ratio := 1.0
        title --> "Field at k=$k"
        fillcolor --> :pu_or

        r = mean_radius(model.particles)/2
        x(t) = r * cos(t) + model.listener_positions[1, 1]
        y(t) = r * sin(t) + model.listener_positions[2, 1]

        (x, y, -π/3, π/3)
    end

end

"Plot the response across all wavenumbers"
@recipe function plot(model::FrequencyModel)
    label --> ["real" "imaginary"]
    xlabel --> "Wavenumber (k)"
    ylabel --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(model.particles[1].r,2)) contained in a $(lowercase(name(model.shape)))\n with volfrac=$(signif(calculate_volfrac(model),2)) measured at ($(model.listener_positions[1,1]), $(model.listener_positions[2,1]))"

    (model.k_arr, [real(model.response) imag(model.response)])
end

"Plot the response across time"
@recipe function plot(model::TimeModel)
    label --> ["real" "imaginary"]
    xlabel --> "Time (t)"
    ylabel --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(model.frequency_model.particles[1].r,2)) contained in a $(lowercase(name(model.frequency_model.shape)))\n with volfrac=$(signif(calculate_volfrac(frequency_model.model),2)) measured at ($(model.frequency_model.listener_positions[1,1]), $(model.frequency_model.listener_positions[2,1]))"

    (model.time_arr, [real(model.response) imag(model.response)])
end

