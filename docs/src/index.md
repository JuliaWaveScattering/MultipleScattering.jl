# MultipleScattering.jl

*A Julia library for simulating, processing, and plotting multiple scattering of waves.*

MultipleScattering.jl calculates how time-harmonic waves are scattered by many particles. It uses the multipole method, also called the T-matrix method, to solve the Helmholtz equation. This method is particularly efficient for particles in an infinite domain. Each particle is represented by its T-matrix, so particles can have any shape and properties for which a T-matrix is known.

The package is set up for different spatial dimensions and different types of waves which satisfy Helmholtz equations, such as acoustics, electromagnetism, and elasticity. Acoustics in 2D and 3D is implemented here. For details on the maths see [Martin (1995)](https://pdfs.semanticscholar.org/8bd3/38ec62affc5c89592a9d6d13f1ee6a7d7e53.pdf), [Gower et al. (2017)](https://arxiv.org/abs/1712.05427), and the notes in [Theory](@ref).

## Features

- Acoustic waves in two and three spatial dimensions.
- Particles: penetrable, sound-hard and sound-soft circles and spheres, concentric capsules, and 2D Helmholtz resonators. New particles only need a T-matrix, see [New particles](@ref).
- Sources (incident waves): plane waves, point sources, beams, and any linear combination of these. See [RegularSources](@ref).
- Results in the frequency domain, and in the time domain by using Fourier transforms with smooth impulses. See [Time response](@ref).
- Random placement of particles in a region, for a given volume fraction, and statistics over many random configurations.
- Plot recipes for particles, shapes, sources, fields, and results. See [Plotting](@ref).

## Installation

Install Julia v1.10 or later, then run

```julia
using Pkg
Pkg.add("MultipleScattering")
```

## Quick example

Calculate the field scattered by a single circular particle when it is hit by a plane wave:

```jldoctest index
using MultipleScattering

host_medium = Acoustic(2; ρ = 1.0, c = 1.0)  # 2D acoustics, density 1.0 and sound speed 1.0
particle = Particle(Acoustic(2; ρ = 0.5, c = 0.6), Sphere([0.0, 0.0], 1.0))
source = plane_source(host_medium; direction = [1.0, 0.0])

# the total field at the point x = [3.0, 0.0] for the angular frequency ω = 2.0
result = run([particle], source, [3.0, 0.0], 2.0)
round.(field(result), digits = 6)

# output

1×1 Matrix{ComplexF64}:
 0.391303 + 1.054422im
```

If you have the package [Plots](http://docs.juliaplots.org) installed, then `plot(FrequencySimulation([particle], source), 2.0)` plots the whole field around the particle.

## Where to start

- The [Introduction](@ref) walks through a complete simulation, from the medium to plotting results in time.
- [Simulations and results](@ref) explains what `run` calculates, how to choose `basis_order`, and how to use the results.
- The manual pages on [RegularSources](@ref), [Shapes](@ref), and [Particles](@ref) describe the building blocks in more detail.
- The [Examples](example/README.md) show complete scripts, many with plots.
- The [Library](@ref base_physical_property) lists every type and function.

## Related packages

These packages build on MultipleScattering.jl:

- [EffectiveWaves.jl](https://github.com/JuliaWaveScattering/EffectiveWaves.jl): ensemble averaged (effective) waves in materials filled with random particles.
- [ElasticWaves.jl](https://github.com/JuliaWaveScattering/ElasticWaves.jl): elastic waves, including scattering from elastic particles.
- [ParticleCorrelations.jl](https://github.com/arturgower/ParticleCorrelations.jl): pair correlations and structure factors of particle configurations.
- [MethodOfFundamentalSolutions.jl](https://github.com/arturgower/MethodOfFundamentalSolutions.jl): the method of fundamental solutions for wave and elastostatic problems.

## Contents

```@contents
Pages = [
    "manual/intro.md",
    "manual/simulation.md",
    "manual/source.md",
    "manual/shapes.md",
    "manual/particles.md",
    "manual/time_response.md",
    "manual/plot.md",
    "manual/new_types.md",
]
Depth = 1
```
