"Build a 'field model' with lots of listeners using the same domain as model for a single wavenumber. This 'field model' can then be used to plot the whole field for this wavenumber."
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

"Plots the field from a model with a grid structure of listeners"
function plot_field_model{T}(field_model::FrequencyModel{T}, k::T; res=10, xres=res, yres=res, resp_fnc=real, zero_centred=true)

    # Build the listeners or pixels
    bounds = bounding_box(field_model.shape)
    box_size = bounds.topright - bounds.bottomleft

    # For this we sample at the centre of each pixel
    x_pixels = linspace(bounds.bottomleft[1], bounds.topright[1], xres+1)
    y_pixels = linspace(bounds.bottomleft[2], bounds.topright[2], yres+1)

    # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
    response_mat = transpose(reshape(field_model.response, (xres+1, yres+1)))

    fillcolor = :auto
    if zero_centred
        max_value_to_plot = maximum(resp_fnc(response_mat))
        min_value_to_plot = minimum(resp_fnc(response_mat))
        zero_point = -min_value_to_plot / (max_value_to_plot - min_value_to_plot)
        if (zero_point <= 1) && (zero_point >= 0)
            # If we can, produce a nice zero centred colour scheme, this should work but it doesn't
            WaveyColours = ColorGradient(
                [
                    RGBA{Float64}(1, 0, 0, 1),
                    RGBA{Float64}(1, 1, 1, 1),
                    RGBA{Float64}(0, 0, 1, 1)
                ],
                [0.0, zero_point, 1.0]
            )
            fillcolor = cgrad(WaveyColours)
        end
    end

    # Finally, draw the field
    plot(x_pixels, y_pixels, resp_fnc(response_mat), linetype=:contour, fill=true, fillcolor=fillcolor)

end

"Plot the field for a specific wavenumber from a model you have already set up."
function plot_field{T}(model::FrequencyModel{T}, k::T; res=10, xres=res, yres=res, resp_fnc=real, zero_centred=true, domain_flag=true, particles_flag=domain_flag, shape_flag=domain_flag, listeners_flag=domain_flag)
    field_model = build_field_model(model, k; xres=xres, yres=yres)

    plot_field_model(field_model, k; xres=xres, yres=yres, resp_fnc=resp_fnc, zero_centred=zero_centred)

    # Finally, add particles and listener position from original model
    plot_domain(model; newplot=false, particles_flag=particles_flag, shape_flag=shape_flag, listeners_flag=listeners_flag)

end
