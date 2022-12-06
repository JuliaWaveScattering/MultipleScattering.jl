__precompile__()

module MultipleScattering

## Symmetries
export Symmetry, AbstractSymmetry, AbstractPlanarSymmetry, AbstractAzimuthalSymmetry
export WithoutSymmetry, PlanarSymmetry, PlanarAzimuthalSymmetry, AzimuthalSymmetry, RadialSymmetry, TranslationSymmetry

## Shapes
export Shape, Circle, Sphere, Rectangle, Box, EmptyShape, Halfspace, Plate, TimeOfFlightPlaneWaveToPoint, TimeOfFlightPointWaveToPoint

export outer_radius, volume, name, iscongruent, (≅), congruent, in, issubset, origin, shape, (==), isequal, show
export boundary_functions, boundary_points, boundary_data, bounding_box, corners
export points_in_shape, bottomleft, topright

## Physical mediums
export PhysicalMedium, outgoing_basis_function, regular_basis_function, outgoing_radial_basis, regular_radial_basis, outgoing_translation_matrix, regular_translation_matrix, estimate_regular_basisorder, estimate_outgoing_basisorder, basisorder_to_basislength, basis_coefficients, basislength_to_basisorder, internal_field, check_material
export spatial_dimension, field_dimension

## Electromagnetic
export Electromagnetic

## Acoustics
export Acoustic, AcousticCapsule, sound_hard, hard, rigid, zero_neumann, sound_soft, soft, zero_dirichlet, pressure_release, impedance

## Particles
export AbstractParticle, Particle, CapsuleParticle, AbstractParticles

## RegularSources
export AbstractSource, RegularSource, source_expand, regular_spherical_coefficients, self_test, constant_source, point_source, plane_source, regular_spherical_source, (*), (+)
export PlaneSource

## Main simulation and results
export Simulation, run, TimeSimulation, forcing, field
export SimulationResult, FrequencySimulation, FrequencySimulationResult, size


## Impulses for time response
export ContinuousImpulse, TimeDiracImpulse, FreqDiracImpulse, GaussianImpulse
export DiscreteImpulse, continuous_to_discrete_impulse, DiscreteTimeDiracImpulse, DiscreteGaussianImpulse
export TimeSimulationResult, frequency_to_time, time_to_frequency
export ω_to_t, t_to_ω, firstnonzero

export random_particles, statistical_moments

export t_matrix, get_t_matrices
export scattering_matrix

import Printf: @printf

using StaticArrays: SVector
using OffsetArrays: OffsetArray

using SpecialFunctions: besselj, hankelh1
using WignerSymbols, GSL
using Random, LinearAlgebra, RecipesBase, Statistics
using ProgressMeter


# Generic machinery common to all physical models
include("types.jl")
include("shapes/shapes.jl")

include("physics/special_functions.jl") # Special functions missing from Base library
include("physics/physical_medium.jl")

include("particle.jl")
include("source.jl")
include("result.jl")
include("simulation.jl")
include("impulse.jl")
include("time_simulation.jl")
include("random/random.jl")
include("t_matrix.jl")
include("scattering_matrix.jl")

# Specific physical models
include("physics/diffbessel.jl")
include("physics/acoustics/export.jl")
include("physics/electromagnetism.jl")

#Plotting and graphics
include("plot/plot.jl")

# Precompile hints
include("precompile.jl")

end # module
