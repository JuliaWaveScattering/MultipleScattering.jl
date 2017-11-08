using MultipleScattering
using MultipleScattering.Plot

Plots.plot()

rectangle = Rectangle([0.0,0.0],[2.0,3.0])
plot_shape(rectangle)

circle = Circle(3.0,[-1.0,2.0])
plot_shape(circle)

timeofflight = TimeOfFlight([-1.0,0.0],3.0)
plot_shape(timeofflight)
