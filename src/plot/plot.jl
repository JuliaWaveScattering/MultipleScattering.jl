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

"Plot the field for a particular wavenumber"
@recipe function plot(sim::FrequencySimulation, ω::Number; res=10, xres=res, yres=res,
                         field_apply=real, build_field=true, bounds = :auto,
                         drawparticles=false)

    @series begin
        # find a box which covers everything
        if bounds == :auto bounds = bounding_rectangle(sim.particles) end

        if build_field
          field_sim = run(sim, bounds, [ω]; xres=xres, yres=yres)
        else
          field_sim = sim
        end

        xy_mat = reshape(field_sim.x, (xres+1, yres+1))

        x_pixels = [x[1] for x in xy_mat[:,1]]
        y_pixels = [x[2] for x in xy_mat[1,:]]

        # Turn the responses (a big long vector) into a matrix, so that the heatmap will understand us
        response_mat = transpose(reshape(field(field_sim), (xres+1, yres+1)))
        linetype --> :contour
        fill --> true
        aspect_ratio := 1.0
        fillcolor --> :pu_or
        title --> "Field at ω=$ω"

        (x_pixels, y_pixels, field_apply.(response_mat))
    end
    if drawparticles # written in strange way due to odd behaviour of @series
      particles = filter(p -> inside(bounds, p), sim.particles)
      for i=1:length(particles) @series particles[i] end
    end
    # if drawsource
    #   @series begin
    #       line --> 0
    #       fill --> (0, :lightgreen)
    #       legend --> false
    #       grid --> false
    #       colorbar --> true
    #       aspect_ratio := 1.0
    #
    #       r = mean_radius(sim.particles)/2
    #       x(t) = r * cos(t) + sim.listener_positions[1, 1]
    #       y(t) = r * sin(t) + sim.listener_positions[2, 1]
    #
    #       (x, y, -2π/3, 2π/3)
    #   end
    # end
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
