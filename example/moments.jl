# Example shows how to egnerate a whole batch of responses with the same label
# (volume fraction, radius and shape) but different realisations (actual
# positions of the particles). We will then extract stochastic information from
# them
using MultipleScattering

volfrac = 0.01
radius = 1.0
k_arr = collect(linspace(0.01,1.0,100))

# Holder for our models
models = Vector{FrequencyModel{Float64}}(10)
for i=1:10
    models[i] = FrequencyModel(volfrac,radius,k_arr)
end

moments = Moments(models)