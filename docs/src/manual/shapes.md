# Shapes

```@meta
DocTestSetup = quote
    using MultipleScattering
end
```
Shape is an abstract type which represents the shape of particles, and also the domain to place particles. See [Shape](@ref base_shape) for a list of relevant types and functions.


## Existing shapes
The package provides three basic shapes. You can plot them using:
```jldoctest intro; output = false
rectangle = Rectangle([0.0,-1.0],[1.0,2.0])
circle = Circle([-1.0,0.0],1.0)
timeofflight = TimeOfFlight([1.0,0.0],3.0)
# output
TimeOfFlight{Float64}([1.0, 0.0], 3.0)
```
```julia
using Plots;
plot(rectangle, linecolor = :red)
plot!(circle, linecolor = :green)
plot!(timeofflight, linecolor = :blue)
```
![Plot the three shapes](../media/shapes.png)

The `Rectangle` and `TimeOfFlight` are usually region where particles are placed. Time of flight is a shape which contains shapes from a half space which take at most `t` time to reach from the listener. The `Circle` is also used to define circular particles.
