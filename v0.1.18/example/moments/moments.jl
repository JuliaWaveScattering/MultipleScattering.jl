# Example shows how to egnerate a whole batch of responses with the same label
# (volume fraction, radius and shape) but different realisations (actual
# positions of the particles). We will then extract stochastic information from
# them
using MultipleScattering

function moments_example()
    volfrac = 0.01
    radius = 1.0
    num_particles = 10
    k_arr = collect(LinRange(0.01,1.0,100))

    # region to place the particles
    shape = Rectangle(volfrac, radius, num_particles)
    # Holder for our simulations
    simulations = Vector{FrequencySimulation{Float64}}(undef,10)
    for i=1:10
        simulations[i] = FrequencySimulation(volfrac,radius,k_arr;seed=[0x7f5def91, 0x82255da3, 0xc5e461c7, UInt32(i)])
    end

    moments = StatisticalMoments(simulations)
    return moments
end
