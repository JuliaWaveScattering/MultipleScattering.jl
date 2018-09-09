__precompile__()

module MultipleScattering

export Shape, Circle, Rectangle, outer_radius, volume, name, iscongruent,
       congruent, in, issubset, origin, shape, Sphere, TimeOfFlight,
       TimeOfFlightFromPoint, (==), isequal
export boundary_functions, boundary_points, boundary_data, bounding_rectangle
export bottomleft, topright
export PhysicalProperties, Acoustic, Electromagnetic, AcousticCapsule,
       basis_function, basis_coefficients, internal_field, dim, field_dim,
       sound_hard, hard, rigid, zero_neumann, sound_soft, soft, zero_dirichlet,
       pressure_release, impedance
export AbstractParticle, Particle, CapsuleParticle, AbstractParticles
export Source, besselj_field, self_test, point_source, plane_source,
       (*), (+)
export Simulation, run, TimeSimulation, forcing, field
export SimulationResult, FrequencySimulation, FrequencySimulationResult, size
export ContinuousImpulse, TimeDiracImpulse, FreqDiracImpulse, GaussianImpulse
export DiscreteImpulse, DiscreteTimeDiracImpulse, DiscreteGaussianImpulse
export TimeSimulationResult, frequency_to_time, time_to_frequency
export ω_to_t, t_to_ω, firstnonzero
export random_particles, calculate_moments
export t_matrix, get_t_matrices
export scattering_matrix

import SpecialFunctions: besselj, hankelh1

import StaticArrays: SVector

import OffsetArrays: OffsetArray

using RecipesBase
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
