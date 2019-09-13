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

## New shape
If you are feeling very adventurous, you can define your own shape. As an example we show how to define a rectangle.
```julia
struct MyRectangle{T} <: MultipleScattering.Shape{T,2}
    origin::Vector{T}
    width::T
    height::T
end
```

To be able to place particles in `MyRectangle` you need to define the characteristics below. To define a new particle using `MyRectangle` you would need to also define
```julia
MultipleScattering.volume(shape::MyRectangle) = shape.width * shape.height

MultipleScattering.name(shape::MyRectangle) = "MyRectangle"

# every shape needs a bounding rectangle, this is where particles are first placed.
MultipleScattering.bounding_rectangle(shape::MyRectangle) =  MultipleScattering.Rectangle(shape.origin, shape.width, shape.height)

import Base.in
function in(x::AbstractVector, r::MyRectangle)::Bool
    all(abs.(x .- r.origin) .<= [r.width, r.height])
end

import Base.issubset
function issubset(circle::MultipleScattering.Circle{T}, r::MyRectangle{T}) where {T}
    all((MultipleScattering.origin(circle) .- circle.radius) .>= r.origin - [r.width/2.0, r.height/2.0]) &&
    all((MultipleScattering.origin(circle) .+ circle.radius) .<= r.origin + [r.width/2.0, r.height/2.0])
end
```
With these definitions you can now fill this shape with circular particles, as shown in [Particles](@ref).
