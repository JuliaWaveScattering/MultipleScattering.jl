# New particles

If you are feeling very adventurous, you can define a new type of particle.

A particle is defined by its shape and physical properties. The simplest example being a circular acoustics particle: `Particle(Acoustic(1.0,1.0,2),Circle(1.0))`. To define a new particle you need to either use a new shape or a new physical medium.

## New shape
Suppose we want to define a new shape, which is very much like a rectangle.

```julia
using MultipleScattering

struct MyRectangle{T} <: Shape{2}
    origin::Vector{T}
    width::T
    height::T
end
```

To be able to place particles in `MyRectangle` you need to extend the functions defined in the package you need to explicitly import them:
```julia
import MultipleScattering: volume, name, bounding_box

volume(shape::MyRectangle) = shape.width * shape.height
name(shape::MyRectangle) = "MyRectangle"

# every shape needs a bounding rectangle, this is where particles are first placed.
bounding_box(shape::MyRectangle) = Box(shape.origin, [shape.width, shape.height])

import Base.in, Base.issubset
function in(x::AbstractVector, r::MyRectangle)::Bool
    all(abs.(x .- r.origin) .<= [r.width, r.height] ./ 2)
end

function issubset(circle::Sphere{T,2}, r::MyRectangle{T}) where {T}
    all((origin(circle) .- circle.radius) .>= r.origin - [r.width/2.0, r.height/2.0]) &&
    all((origin(circle) .+ circle.radius) .<= r.origin + [r.width/2.0, r.height/2.0])
end
```
!!! note
    Using other packages with MultipleScattering may result in a naming conflict for functions like `volume`. In this case you will have to explicitly call `MultipleScattering.volume`.

With these definitions you can now fill this shape with circular particles, as shown in [Placing particles in a region](@ref).

To define a new particle from `MyRectangle` you need to extend the functions [`outer_radius`](@ref), `Base.(==)`, [`iscongruent`](@ref), then supposing you will use a predefined `PhysicalMedium` (such as `Acoustic`), you then need to extend [`t_matrix`](@ref).

## New T-matrix

The scattering from a particle is completely determined by its T-matrix. To simulate a new type of particle, you need to extend [`t_matrix`](@ref) for that particle and the background medium. The T-matrix maps the coefficients of the wave incident on the particle, in the regular basis, to the coefficients of the scattered wave, in the outgoing basis. See [What `run` calculates](@ref run_maths) for the conventions, and the [Notes on the T-matrix](../maths/multiplescattering.pdf).

For example, in 2D a circular particle where the field is zero on the boundary (a sound-soft particle), with radius $a$, has the diagonal T-matrix with entries $-\mathrm J_n(k a) / \mathrm H^{(1)}_n(k a)$, for $n = -N, \ldots, N$, where $N$ is the basis order.

## New physical medium

A physical medium is a subtype of [`PhysicalMedium`](@ref)`{Dim,FieldDim}`, where `Dim` is the spatial dimension and `FieldDim` is the dimension of the field, for example `FieldDim = 1` for a scalar field.

For a scalar wave which satisfies a Helmholtz equation, the medium only needs a field `c`, the wave speed. Then the basis functions and translation matrices of the package are used. You need to define the T-matrix of the particles, and a source. For example, the code below defines a new medium and sound-soft circular particles in that medium:

```julia
using MultipleScattering, LinearAlgebra, SpecialFunctions
import MultipleScattering: t_matrix

struct MyMedium{T,Dim} <: PhysicalMedium{Dim,1}
    c::Complex{T} # the wave speed
end

function t_matrix(p::Particle{2,MyMedium{T,2},Sphere{T,2}}, medium::MyMedium{T,2}, ω::T, basis_order::Integer) where T
    ka = ω / medium.c * outer_radius(p)
    return Diagonal([- besselj(n, ka) / hankelh1(n, ka) for n = -basis_order:basis_order])
end

medium = MyMedium{Float64,2}(1.0)

# a source given by its coefficients in the regular basis, here the wave J_0(k r)
source = regular_spherical_source(medium, [1.0 + 0.0im])

particles = [Particle(medium, Sphere([2.0, 0.0], 0.5)), Particle(medium, Sphere([-1.0, 1.0], 0.5))]

result = run(particles, source, [[4.0, 0.0], [0.0, 0.0]], 1.0)
```

This gives the same result as using acoustic particles with the medium `sound_soft(2)` in the medium `Acoustic(2; ρ = 1.0, c = 1.0)`. To also calculate the field inside the particles, you need to extend [`internal_field`](@ref).

For waves with more than one component, such as elastic waves, you also need to extend the functions which define the basis:

- [`basisorder_to_basislength`](@ref) and [`basislength_to_basisorder`](@ref), the number of basis functions for each basis order;
- [`regular_basis_function`](@ref) and [`outgoing_basis_function`](@ref);
- [`regular_translation_matrix`](@ref) and [`outgoing_translation_matrix`](@ref).

See the package [ElasticWaves.jl](https://github.com/JuliaWaveScattering/ElasticWaves.jl) for an example.
