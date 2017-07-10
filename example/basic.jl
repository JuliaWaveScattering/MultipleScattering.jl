# If it isn't installed, clone it from github
try using MultipleScattering
catch
    Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
end

using MultipleScattering

volfrac = 0.01
radius = 1.0
k_arr = collect(linspace(0.01,1.0,100))
model = FrequencyModel(volfrac,radius,k_arr)

using MultipleScattering.Plot

plot_model(model)

using PyPlot
PyPlot.figure(2)

plot_field(model,0.5)