using MultipleScattering

k_arr = collect(linspace(0.1,1.0,10))

# You can also pick your own shape, an generate random particles inside it 
# with a certain radius ands volume fraction
radius = 0.5
volfrac = 0.2

circle = Circle(5.0,[0.0,0.0])
circle_model = FrequencySimulation(volfrac,radius,k_arr;shape=circle)

using Plots
plot(circle_model,0.5)
