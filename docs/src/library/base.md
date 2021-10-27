# Base

```@meta
CurrentModule = MultipleScattering
```

```@contents
Pages = ["base.md"]
```

## [Physical mediums](@id base_physical_property)

The types and functions shared between all different types of physical mediums.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/physics/physical_medium.jl", "src/physics/special_functions.jl"]
```

## [Shapes and domains](@id base_shapes)

Shapes are used to define the shape of particles, and to define containers (or configurations) where all particles are placed. It is also used in plotting.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/shapes/shapes.jl", "src/shapes/sphere.jl", "src/shapes/halfspace.jl", "src/shapes/plate.jl", "src/shapes/box.jl", "src/shapes/empty_shape.jl", "src/shapes/time_of_flight.jl"]
```

## [Particle](@id base_particle)

Particle types and functions.

```@docs
AbstractParticle
Particle
CapsuleParticle
iscongruent(::AbstractParticle,::AbstractParticle)
```

## [RegularSource](@id source_base)

RegularSource types and functions for all physical mediums.

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/source.jl"]
```

## [Impulse](@id impulse_base)

For convert to the time domain.

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
