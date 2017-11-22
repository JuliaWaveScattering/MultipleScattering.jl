using MultipleScattering
using Plots; pyplot()

function hankel_order_convergence(m=[0,1,2,3,4,5,6,7,8,9,10], volfrac = 0.1,
    radius = 1.0, maxtime = 40.0, k_arr=collect(linspace(0.01,1.0,100)) )

    listener_position = [-10.0,0.0]
    shape = TimeOfFlight(listener_position,maxtime)

    seed = MersenneTwister(1).seed
    particles = random_particles(volfrac, radius, shape; seed = seed)

    models = Vector{FrequencyModel{Float64}}(length(m))

    for i = eachindex(m)
        models[i] = FrequencyModel(particles, k_arr; seed=seed, hankel_order=m[i])
    end
    
    return models
end

function plot_hankel_order_convergence(models)
    responses = Vector{Vector{Complex{Float64}}}(length(models))
    m = Vector{Int64}(length(models))

    labels = Matrix{String}(1,0)
    for i = eachindex(models)
        responses[i] = reshape(models[i].response, size(models[i].response,1))
        m[i] = models[i].hankel_order
        labels = [labels "m = $(m[i])"]
    end

    error = [r .- responses[end] for r in responses[1:end-1]]
    integrated_error = norm.(error).*map(m->((m.k_arr[end]-m.k_arr[1])/length(m.k_arr)),models[1:end-1])

    colors = reshape(linspace(RGB(0.6,1,0.6),RGB(0,0.4,0),length(m)),1,length(m))
    realcolors = reshape(linspace(RGB(0.6,0.6,1),RGB(0,0,0.4),length(m)),1,length(m))
    imagcolors = reshape(linspace(RGB(1,0.6,0.6),RGB(0.4,0,0),length(m)),1,length(m))

    absvec(v) = abs.(v)
    plot(
        plot(models[end],0.5),
        plot(models[1].k_arr, [real(responses),imag(responses)],
             labels=[("real ".*labels) ("imag ".*labels)],
             xlab="Wavenumber (k)", ylab="Response", linecolor=[realcolors imagcolors]
        ),
        plot(models[1].k_arr, absvec.(v),
             yscale=:log10, labels=labels, linecolor=colors, 
             xlab="Wavenumber (k)", ylab="Absolute error",
        ),
        plot(m[1:end-1], integrated_error,
             yscale=:log10, legend=false,
             xlab="Hankel order", ylab="\$L^2\$ integrated error",
        )
    )

end

models = hankel_order_convergence()
plot_hankel_order_convergence(models)
