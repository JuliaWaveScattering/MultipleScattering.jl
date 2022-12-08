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
