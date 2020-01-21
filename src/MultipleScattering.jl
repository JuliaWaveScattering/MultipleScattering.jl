__precompile__()

module MultipleScattering

## Shapes
export Shape, Circle, Rectangle, EmptyShape, Halfspace, TimeOfFlight, TimeOfFlightFromPoint

export outer_radius, volume, name, iscongruent, (≅), congruent, in, issubset, origin, shape, Sphere, (==), isequal, show
export boundary_functions, boundary_points, boundary_data, bounding_rectangle
export bottomleft, topright

## Physical mediums
export PhysicalMedium, Acoustic, Electromagnetic, AcousticCapsule, outgoing_basis_function, regular_basis_function, basis_coefficients, internal_field, sound_hard, hard, rigid, zero_neumann, sound_soft, soft, zero_dirichlet, pressure_release, impedance
export spatial_dimension, field_dimension

## Particles
export AbstractParticle, Particle, CapsuleParticle, AbstractParticles

## Sources
export AbstractSource, Source, source_expand, self_test, constant_source, point_source, plane_source,
       (*), (+)

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

import SpecialFunctions: besselj, hankelh1
import Printf: @printf

import StaticArrays: SVector
import OffsetArrays: OffsetArray

using Random, LinearAlgebra, RecipesBase, Statistics
using ProgressMeter

# Generic machinery common to all physical models
include("shapes/shapes.jl")
include("physics/physical_properties.jl")
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
include("physics/acoustics/acoustics.jl")
include("physics/electromagnetism.jl")

#Plotting and graphics
include("plot/plot.jl")

# Precompile hints
include("precompile.jl")

end # module
