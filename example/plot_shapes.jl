using MultipleScattering
using Plots

rectangle = Rectangle([0.0,0.0],[2.0,3.0])
plot(rectangle)

circle = Circle(3.0,[-1.0,2.0])
plot!(circle)

timeofflight = TimeOfFlight([-1.0,0.0],3.0)
plot!(timeofflight)
