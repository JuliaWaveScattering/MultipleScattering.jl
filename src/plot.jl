include("plot_domain.jl")
include("plot_moments.jl")

"""
Build a 'field simulation' with lots of listeners using the same domain as simulation
you pass in. This 'field simulation' can then be used to plot the whole field for
this wavenumber.
"""
function build_field_simulation{T}(simulation::FrequencySimulation{T}, bounds::Rectangle{T},
                              k_arr::Vector{T}=simulation.k_arr; res=10,xres=res,yres=res)
    # Create deep copy of simulation so that we can add lots of new listener positions and rerun the simulation
    field_simulation = deepcopy(simulation)

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

    field_simulation.listener_positions = listener_positions
    generate_responses!(field_simulation,k_arr)

    return field_simulation
end

"Plot the field for a particular wavenumber"
@recipe function plot{T}(simulation::FrequencySimulation{T},k::T;res=10, xres=res, yres=res,
                         resp_fnc=real, build_field=true,
                         drawparticles=true, drawshape=false, drawlisteners=build_field)

    @series begin
        # find a box which covers everything
        shape_bounds = bounding_box(simulation.shape)
        listeners_as_particles = map(
            l -> Particle(simulation.listener_positions[:,l],mean_radius(simulation)/2),
            1:size(simulation.listener_positions,2)
        )
        particle_bounds = bounding_box([simulation.particles; listeners_as_particles])
        bounds = bounding_box(shape_bounds, particle_bounds)
        if build_field
          field_simulation = build_field_simulation(simulation, bounds, [k]; xres=xres, yres=yres)
        else
          field_simulation = simulation
        end

        # For this we sample at the centre of each pixel
        x_pixels = linspace(bounds.bottomleft[1], bounds.topright[1], xres+1)
        y_pixels = linspace(bounds.bottomleft[2], bounds.topright[2], yres+1)

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field_simulation.response, (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        fillcolor --> :pu_or
        title --> "Field at k=$k"

        (x_pixels, y_pixels, resp_fnc.(response_mat))
    end
    if drawshape
      @series begin
          simulation.shape
      end
    end
    if drawparticles
      for i=1:length(simulation.particles) @series simulation.particles[i] end
    end  
    if drawlisteners
      @series begin
          line --> 0
          fill --> (0, :lightgreen)
          legend --> false
          grid --> false
          colorbar --> true
          aspect_ratio := 1.0

          r = mean_radius(simulation.particles)/2
          x(t) = r * cos(t) + simulation.listener_positions[1, 1]
          y(t) = r * sin(t) + simulation.listener_positions[2, 1]

          (x, y, -2π/3, 2π/3)
      end
    end
end

"Plot the response across all wavenumbers"
@recipe function plot(simulation::FrequencySimulation)
    label --> ["real" "imaginary"]
    xlabel --> "Wavenumber (k)"
    ylabel --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(simulation.particles[1].r,2)) contained in a $(lowercase(name(simulation.shape)))\n with volfrac=$(signif(calculate_volfrac(simulation),2)) measured at ($(simulation.listener_positions[1,1]), $(simulation.listener_positions[2,1]))"

    (simulation.k_arr, [real(simulation.response) imag(simulation.response)])
end

"Plot the response across time"
@recipe function plot(simulation::TimeSimulation)
    xlabel --> "Time (t)"
    ylabel --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(simulation.frequency_simulation.particles[1].r,2)) contained in a $(lowercase(name(simulation.frequency_simulation.shape)))\n with volfrac=$(signif(calculate_volfrac(simulation.frequency_simulation),2)) measured at ($(simulation.frequency_simulation.listener_positions[1,1]), $(simulation.frequency_simulation.listener_positions[2,1]))"

    (simulation.time_arr, simulation.response)
end

"Plot the field for a particular an array of time"
@recipe function plot{T}(TimeSimulation::TimeSimulation{T}, t::Union{T,Vector{T}};
                        res=10, xres=res, yres=res, resp_fnc=real,
                        drawshape = false, drawlisteners = false)
    simulation = TimeSimulation.frequency_simulation

    if isa(t,T) t_arr = [t]
    else        t_arr = t  end

    @series begin
        # find a box which covers everything
        shape_bounds = bounding_box(simulation.shape)
        listeners_as_particles = map(
            l -> Particle(simulation.listener_positions[:,l],mean_radius(simulation)/2),
            1:size(simulation.listener_positions,2)
        )
        particle_bounds = bounding_box([simulation.particles; listeners_as_particles])
        bounds = bounding_box(shape_bounds, particle_bounds)

        field_simulation = build_field_simulation(simulation, bounds; xres=xres, yres=yres)
        field_TimeSimulation = deepcopy(TimeSimulation) # to use all the same options/fields as TimeSimulation
        field_TimeSimulation.frequency_simulation = field_simulation
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
          simulation.shape
      end
    end
    for i=1:length(simulation.particles) @series simulation.particles[i] end
    if drawlisteners
      @series begin
          line --> 0
          fill --> (0, :lightgreen)
          grid --> false
          colorbar --> true
          aspect_ratio --> 1.0
          legend --> false

          r = mean_radius(simulation.particles)/2
          x(t) = r * cos(t) + simulation.listener_positions[1, 1]
          y(t) = r * sin(t) + simulation.listener_positions[2, 1]

          (x, y, -2π/3, 2π/3)
      end
    end

    # clims --> (-0.65,0.65) # so that white is always = 0, these values are suitable when using default gaussian impulse
end
