# Base

```@meta
CurrentModule = MultipleScattering
```

```@contents
Pages = ["base.md"]
```

## [Shapes](@id base_shapes)

Shape types and functions.

```@docs
Shape
origin
iscongruent(::Shape,::Shape)
congruent
bounding_rectangle
boundary_functions
name
outer_radius
volume
Circle
Rectangle
TimeOfFlightFromPoint
TimeOfFlight
Sphere
```

## Physical properties

Physical properties types and functions.

```@docs
PhysicalProperties
field_dim
dim
outgoing_basis_function
regular_basis_function
internal_field
boundary_data
```

## Particles

Particle types and functions.

```@docs
AbstractParticle
Particle
CapsuleParticle
iscongruent(::AbstractParticle,::AbstractParticle)
```

## Source

Source types and functions.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/source.jl"]
```

## Impulse

Impulse types and functions.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/impulse.jl"]
```


## Simulation

Simulation types and functions.

```@docs
FrequencySimulation
run(::FrequencySimulation)
run(::FrequencySimulation, ::Rectangle, ::AbstractVector)
FrequencySimulationResult
forcing
basis_coefficients
field
scattering_matrix
t_matrix
get_t_matrices
```
