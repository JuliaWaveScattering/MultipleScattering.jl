using MultipleScattering
using Plots; pyplot()

function hankel_order_convergence(m=[0,1,2,3,4,5,6,7,8,9,10], volfrac = 0.1,
    radius = 1.0, maxtime = 40.0, k_arr=collect(LinRange(0.01,1.0,100)) )

    listener_position = [-10.0,0.0]
    shape = TimeOfFlightPlaneWaveToPoint(listener_position,maxtime)

    seed = MersenneTwister(1).seed
    particles = random_particles(volfrac, radius, shape; seed = seed)

    simulations = Vector{FrequencySimulation{Float64}}(undef,length(m))

    for i = eachindex(m)
        simulations[i] = FrequencySimulation(particles, k_arr; seed=seed, hankel_order=m[i])
    end

    return simulations
end

function plot_hankel_order_convergence(simulations)
    responses = Vector{Vector{Complex{Float64}}}(undef,length(simulations))
    m = Vector{Int64}(undef,length(simulations))

    # labels = Matrix{String}(undef,1,0)
    for i = eachindex(simulations)
        responses[i] = reshape(simulations[i].response, size(simulations[i].response,1))
        m[i] = simulations[i].hankel_order
        labels = [labels "m = $(m[i])"]
    end

    error = [r .- responses[end] for r in responses[1:end-1]]
    integrated_error = norm.(error).*map(m->((m.k_arr[end]-m.k_arr[1])/length(m.k_arr)),simulations[1:end-1])

    colors = reshape(LinRange(RGB(0.6,1,0.6),RGB(0,0.4,0),length(m)),1,length(m))
    realcolors = reshape(LinRange(RGB(0.6,0.6,1),RGB(0,0,0.4),length(m)),1,length(m))
    imagcolors = reshape(LinRange(RGB(1,0.6,0.6),RGB(0.4,0,0),length(m)),1,length(m))

    absvec(v) = abs.(v)
    plot(
        plot(simulations[end],0.5),
        plot(simulations[1].k_arr, [real(responses),imag(responses)],
             labels=[ map(c -> "real "*c,labels)  map(c -> "imag "*c,labels)],
             xguide ="Wavenumber (k)", yguide ="Response", linecolor=[realcolors imagcolors]
        ),
        plot(simulations[1].k_arr, absvec.(error),
             yscale=:log10, labels=labels, linecolor=colors,
             xguide ="Wavenumber (k)", yguide ="Absolute error",
        ),
        plot(m[1:end-1], integrated_error,
             yscale=:log10, legend=false,
             xguide ="Hankel order", yguide ="\$L^2\$ integrated error",
        )
    )

end
