# If it isn't installed, clone it from github
try using MultipleScattering
catch
    Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
end

using MultipleScattering

# Define the fraction of the volume the particles will take up, their radius and
# which wavenumbers (k) to evaulate at
volfrac = 0.01
radius = 1.0
k_arr = collect(linspace(0.01,1.0,100))
model = FrequencyModel(volfrac,radius,k_arr)

using MultipleScattering.Plot
Plots.pyplot()

plot_model(model)

import PyPlot.figure
figure(2)

plot_field(model,0.8;res=100)