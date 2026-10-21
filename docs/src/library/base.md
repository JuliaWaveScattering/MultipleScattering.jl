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
Pages   = ["src/physics/physical_medium.jl", "src/physics/special_functions.jl", "src/physics/electromagnetism.jl"]
```

## [Particle](@id base_particle)

Particle types and functions.

```@docs
AbstractParticle
```

```@autodocs
Modules = [MultipleScattering]
Order   = [:type, :function]
Pages   = ["src/particle.jl"]
```

## [Symmetries](@id base_symmetries)

Symmetries shared between the source and the material or particles. These can be used to simplify calculations.

```@docs
Symmetry
AbstractSymmetry
WithoutSymmetry
AbstractPlanarSymmetry
PlanarSymmetry
PlanarAzimuthalSymmetry
AbstractAzimuthalSymmetry
AzimuthalSymmetry
RadialSymmetry
MultipleScattering.AbstractTranslationSymmetry
TranslationSymmetry
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
overlapping_pairs
SimulationResult
FrequencySimulationResult
TimeSimulationResult
basis_coefficients
field
scattering_matrix
t_matrix
get_t_matrices
```
