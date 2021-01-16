# New particles

If you are feeling very adventurous, you can define a new type of particle.

A particle is define by its shape and physical properties. The simplest example being a circular acoustics particle: `Particle(Acoustic(1.0,1.0,2),Circle(1.0))`. To define a new particle you need to either use a new shape or a new physical medium.

## New shape
Suppose we want to define a new shape, which is very much like a rectangle.

```julia
using MultipleScattering

struct MyRectangle{T} <: Shape{T,2}
    origin::Vector{T}
    width::T
    height::T
end
```

To be able to place particles in `MyRectangle` you need to extend the functions defined in the package you need to explicitly import them:
```julia
import MultipleScattering: volume, name, bounding_rectangle

volume(shape::MyRectangle) = shape.width * shape.height
name(shape::MyRectangle) = "MyRectangle"

# every shape needs a bounding rectangle, this is where particles are first placed.
bounding_rectangle(shape::MyRectangle) = Rectangle(shape.origin, shape.width, shape.height)

import Base.in, Base.issubset
function in(x::AbstractVector, r::MyRectangle)::Bool
    all(abs.(x .- r.origin) .<= [r.width, r.height])
end

function issubset(circle::Sphere{T,2}, r::MyRectangle{T}) where {T}
    all((origin(circle) .- circle.radius) .>= r.origin - [r.width/2.0, r.height/2.0]) &&
    all((origin(circle) .+ circle.radius) .<= r.origin + [r.width/2.0, r.height/2.0])
end
```
!!! note
    Using other packages with MultipleScattering may result in a naming conflict for functions like `volume`. In this case you will have to explicitly call `MultipleScattering.volume`.

With these definitions you can now fill this shape with circular particles, as shown in [Placing particles in a region](@ref).

To define a new particle from `MyRectangle` you need to extend the functions [`outer_radius`](@ref), `Base.(==)`, [`iscongruent`](@ref), then supposing you will use a predefined `PhysicalProperty` (such as `Acoustic`), you then need to extend [`t_matrix`](@ref).
