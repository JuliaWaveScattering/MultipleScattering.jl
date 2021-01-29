"Plot the response across time"
@recipe function plot(simulation::TimeSimulation)
    xguide --> "Time (t)"
    yguide --> "Response"
    grid --> false
    title --> "Response from particles of radius $(signif(simulation.frequency_simulation.particles[1].r,2)) contained in a $(lowercase(name(simulation.frequency_simulation.shape)))\n with volfrac=$(signif(calculate_volfrac(simulation.frequency_simulation),2)) measured at ($(simulation.frequency_simulation.listener_positions[1,1]), $(simulation.frequency_simulation.listener_positions[2,1]))"

    (simulation.time_arr, simulation.response)
end

# "Plot the field for an array of time"
# @recipe function plot{T}(TimeSimulation::TimeSimulation{T}; res=res)
#     map(axes(TimeSimulation.time_arr,1)) do i
#         simulation = deepcopy(TimeSimulation)
#         simulation.response = reshape(TimeSimulation.response[i,:],1,:)
#         plot(simulation,TimeSimulation.time_arr[i]; res=res, build_field=false)
#     end
# end

"Plot the field for one time"
@recipe function plot{T}(TimeSimulation::TimeSimulation{T}, t::Union{T,Vector{T}};
                        res=20, xres=res, yres=res, resp_fnc=real,
                        drawshape = false, build_field=true, drawlisteners = build_field)
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

        if build_field
          field_simulation = run(simulation, bounds; xres=xres, yres=yres)
          field_TimeSimulation = deepcopy(TimeSimulation) # to use all the same options/fields as TimeSimulation
          field_TimeSimulation.frequency_simulation = field_simulation
          generate_responses!(field_TimeSimulation, t_arr)
        else
          field_TimeSimulation = TimeSimulation
        end

        # For this we sample at the centre of each pixel
        x_pixels = LinRange(bounds.bottomleft[1], bounds.topright[1], xres+1)
        y_pixels = LinRange(bounds.bottomleft[2], bounds.topright[2], yres+1)

        # NOTE only plots the first time plot for now...
        response_mat = transpose(reshape(field_TimeSimulation.response[1,:], (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        fillcolor --> :balance
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
