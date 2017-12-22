include("plot_domain.jl")
include("plot_moments.jl")

"""
Build a 'field model' with lots of listeners using the same domain as model
you pass in. This 'field model' can then be used to plot the whole field for
this wavenumber.
"""
function build_field_model{T}(model::FrequencySimulation{T}, bounds::Rectangle{T},
                              k_arr::Vector{T}=model.k_arr; res=10,xres=res,yres=res)
    # Create deep copy of model so that we can add lots of new listener positions and rerun the model
    field_model = deepcopy(model)

    # Build the listeners or pixels
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
    generate_responses!(field_model,k_arr)

    return field_model
end

"Plot the field for a particular wavenumber"
@recipe function plot{T}(model::FrequencySimulation{T},k::T;res=10, xres=res, yres=res,
                         resp_fnc=real, drawshape = false)

    @series begin
        # find a box which covers everything
        shape_bounds = bounding_box(model.shape)
        listeners_as_particles = map(
            l -> Particle(model.listener_positions[:,l],mean_radius(model)/2),
            1:size(model.listener_positions,2)
        )
        particle_bounds = bounding_box([model.particles; listeners_as_particles])
        bounds = bounding_box(shape_bounds, particle_bounds)
        field_model = build_field_model(model, bounds, [k]; xres=xres, yres=yres)

        # For this we sample at the centre of each pixel
        x_pixels = linspace(bounds.bottomleft[1], bounds.topright[1], xres+1)
        y_pixels = linspace(bounds.bottomleft[2], bounds.topright[2], yres+1)

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field_model.response, (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        fillcolor --> :pu_or
        title --> "Field at k=$k"

        (x_pixels, y_pixels, resp_fnc.(response_mat))
    end
    if drawshape
      @series begin
          model.shape
      end
    end
    for i=1:length(model.particles) @series model.particles[i] end

    @series begin
        line --> 0
        fill --> (0, :lightgreen)
        legend --> false
        grid --> false
        colorbar --> true
        aspect_ratio := 1.0

        r = mean_radius(model.particles)/2
        x(t) = r * cos(t) + model.listener_positions[1, 1]
        y(t) = r * sin(t) + model.listener_positions[2, 1]

        (x, y, -2π/3, 2π/3)
    end

end

"Plot the response across all wavenumbers"
@recipe function plot(model::FrequencySimulation)
    label --> ["real" "imaginary"]
    xlabel --> "Wavenumber (k)"
    ylabel --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(model.particles[1].r,2)) contained in a $(lowercase(name(model.shape)))\n with volfrac=$(signif(calculate_volfrac(model),2)) measured at ($(model.listener_positions[1,1]), $(model.listener_positions[2,1]))"

    (model.k_arr, [real(model.response) imag(model.response)])
end

"Plot the response across time"
@recipe function plot(model::TimeSimulation)
    xlabel --> "Time (t)"
    ylabel --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(model.frequency_model.particles[1].r,2)) contained in a $(lowercase(name(model.frequency_model.shape)))\n with volfrac=$(signif(calculate_volfrac(model.frequency_model),2)) measured at ($(model.frequency_model.listener_positions[1,1]), $(model.frequency_model.listener_positions[2,1]))"

    (model.time_arr, model.response)
end

"Plot the field for a particular an array of time"
@recipe function plot{T}(TimeSimulation::TimeSimulation{T}, t_arr;
                        res=10, xres=res, yres=res, resp_fnc=real,
                        drawshape = false, drawlisteners = false)
    model = TimeSimulation.frequency_model

    @series begin
        # find a box which covers everything
        shape_bounds = bounding_box(model.shape)
        listeners_as_particles = map(
            l -> Particle(model.listener_positions[:,l],mean_radius(model)/2),
            1:size(model.listener_positions,2)
        )
        particle_bounds = bounding_box([model.particles; listeners_as_particles])
        bounds = bounding_box(shape_bounds, particle_bounds)

        field_model = build_field_model(model, bounds; xres=xres, yres=yres)
        field_TimeSimulation = deepcopy(TimeSimulation) # to use all the same options/fields as TimeSimulation
        field_TimeSimulation.frequency_model = field_model
        generate_responses!(field_TimeSimulation, t_arr)

        # For this we sample at the centre of each pixel
        x_pixels = linspace(bounds.bottomleft[1], bounds.topright[1], xres+1)
        y_pixels = linspace(bounds.bottomleft[2], bounds.topright[2], yres+1)


        # NOTE only plots the first time plot for now...
        response_mat = transpose(reshape(field_TimeSimulation.response[1,:], (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        fillcolor --> :pu_or
        title --> "Field at time=$(t_arr[1])"

        (x_pixels, y_pixels, resp_fnc.(response_mat))
    end
    if drawshape
      @series begin
          model.shape
      end
    end
    for i=1:length(model.particles) @series model.particles[i] end
    if drawlisteners
      @series begin
          line --> 0
          fill --> (0, :lightgreen)
          grid --> false
          colorbar --> true
          aspect_ratio --> 1.0
          legend --> false

          r = mean_radius(model.particles)/2
          x(t) = r * cos(t) + model.listener_positions[1, 1]
          y(t) = r * sin(t) + model.listener_positions[2, 1]

          (x, y, -2π/3, 2π/3)
      end
    end

    # clims --> (-0.65,0.65) # so that white is always = 0, these values are suitable when using default gaussian impulse
end
