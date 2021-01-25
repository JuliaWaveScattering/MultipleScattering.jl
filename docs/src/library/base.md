# Base

```@meta
CurrentModule = MultipleScattering
```

```@contents
Pages = ["base.md"]
```

## [Shape](@id base_shape)

Shape types and functions.

```@docs
Shape
origin
iscongruent(::Shape,::Shape)
congruent
bounding_box
boundary_functions
name
outer_radius
volume
Sphere
Box
TimeOfFlightFromPoint
TimeOfFlight
```

## [Physical property](@id base_physical_property)

Physical properties types and functions.

```@docs
PhysicalMedium
field_dimension
spatial_dimension
outgoing_basis_function
regular_basis_function
internal_field
boundary_data
```

## [Particle](@id base_particle)

Particle types and functions.

```@docs
AbstractParticle
Particle
CapsuleParticle
iscongruent(::AbstractParticle,::AbstractParticle)
```

## [Source](@id source_base)

Source types and functions.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/source.jl"]
```

## [Impulse](@id impulse_base)

Impulse types and functions.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/impulse.jl", "src/time_simulation.jl"]
```


## [Simulation](@id simulation_base)

Simulation types and functions.

```@docs
FrequencySimulation
run(::FrequencySimulation)
run(::FrequencySimulation, ::Shape, ::AbstractVector)
FrequencySimulationResult
basis_coefficients
field
scattering_matrix
t_matrix
get_t_matrices
```
