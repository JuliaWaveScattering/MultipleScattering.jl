include("plot_domain.jl")
# include("plot_moments.jl")


"Plot just the particles and source"
@recipe function plot(sim::FrequencySimulation; bounds = :auto,
                         drawparticles=true, drawsource=true)

    if bounds == :auto
      bounds = bounding_rectangle(sim.particles)
      # bounds = bounding_rectangle(shape_bounds, particle_bounds)
    end

    # if drawparticles # written in strange way due to odd behaviour of @series
      particles = filter(p -> inside(bounds, p), sim.particles)
      for i=1:length(particles) @series particles[i] end
    # end
end

"""
Build a 'field simulation' with lots of listeners using the same domain as simulation
you pass in. This 'field simulation' can then be used to plot the whole field for
this wavenumber.
"""
function build_field_simulation{T}(simulation::FrequencySimulation{T}, bounds::Rectangle{T},
                              k_arr::Vector{T}=simulation.k_arr; res=20,xres=res,yres=res)
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
@recipe function plot(sim::FrequencySimulation,k::Number;res=10, xres=res, yres=res,
                         resp_fnc=real, build_field=true, bounds = :auto,
                         drawparticles=true, drawshape=false, drawlisteners=build_field)

# k=0.1; res = 10; xres=res; yres=res; bounds = :auto
# drawparticles=true

    @series begin
        # find a box which covers everything
        if bounds == :auto
          particle_bounds = bounding_rectangle(sim.particles)
          # particle_bounds = bounding_rectangle([sim.particles; listeners_as_particles])
          # bounds = bounding_rectangle(shape_bounds, particle_bounds)
        end
        if build_field
          field_sim = build_field_simulation(sim, bounds, [k]; xres=xres, yres=yres)
        else
          field_sim = sim
        end

        # For this we sample at the centre of each pixel
        x_pixels = linspace(bounds.bottomleft[1], bounds.topright[1], xres+1)
        y_pixels = linspace(bounds.bottomleft[2], bounds.topright[2], yres+1)

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field_sim.response, (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        fillcolor --> :pu_or
        title --> "Field at k=$k"

        (x_pixels, y_pixels, resp_fnc.(response_mat))
    end
    if drawshape
      @series begin
          sim.shape
      end
    end
    if drawparticles # written in strange way due to odd behaviour of @series
      particles = filter(p -> inside(bounds, p), sim.particles)
      for i=1:length(particles) @series particles[i] end
    end
    if drawlisteners
      @series begin
          line --> 0
          fill --> (0, :lightgreen)
          legend --> false
          grid --> false
          colorbar --> true
          aspect_ratio := 1.0

          r = mean_radius(sim.particles)/2
          x(t) = r * cos(t) + sim.listener_positions[1, 1]
          y(t) = r * sin(t) + sim.listener_positions[2, 1]

          (x, y, -2π/3, 2π/3)
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
