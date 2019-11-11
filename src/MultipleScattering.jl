__precompile__()

module MultipleScattering

export Shape, Circle, Rectangle, EmptyShape, outer_radius, volume, name, iscongruent, (≅),
       congruent, in, issubset, origin, shape, Sphere, TimeOfFlight,
       TimeOfFlightFromPoint, (==), isequal, show
export boundary_functions, boundary_points, boundary_data, bounding_rectangle
export bottomleft, topright
export PhysicalProperties, Acoustic, Electromagnetic, AcousticCapsule, outgoing_basis_function,
       regular_basis_function, basis_coefficients, internal_field, dim, field_dim,
       sound_hard, hard, rigid, zero_neumann, sound_soft, soft, zero_dirichlet,
       pressure_release, impedance
export AbstractParticle, Particle, CapsuleParticle, AbstractParticles
export Source, source_expand, self_test, constant_source, point_source, plane_source,
       (*), (+)

export Simulation, run, TimeSimulation, forcing, field
export SimulationResult, FrequencySimulation, FrequencySimulationResult, size

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
