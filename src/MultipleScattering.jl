__precompile__()

module MultipleScattering

export Shape, Circle, Rectangle, outer_radius, volume, name, congruent, inside,
       origin, Sphere, TimeOfFlight, TimeOfFlightFromPoint
export boundary_functions, boundary_points, boundary_data, bounding_rectangle
export bottomleft, topright
export PhysicalProperties, Acoustic, Electromagnetic, AcousticCapsule,
       basis_function, basis_coefficients, inner_basis_coefficients, dim, field_dim
export Particle, (==), Particles
export Source, besselj_field, self_test, TwoDimAcousticPointSource, TwoDimAcousticPlanarSource,
       (*), (+)
export Simulation, run, TimeSimulation, forcing, field
export SimulationResult, FrequencySimulation, FrequencySimulationResult
export TimeSimulationResult, frequency_to_time, time_to_frequency
export ω_to_t, t_to_ω, firstnonzero, get_gaussian_freq_impulse, get_gaussian_time_impulse, delta_freq_impulse
export SimulationDistribution, FrequencySimulationDistribution, sample
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
include("time_simulation.jl")
include("random/random.jl")
include("t_matrix.jl")
include("scattering_matrix.jl")

# Specific physical models
include("physics/diffbessel.jl")
include("physics/acoustics.jl")
include("physics/electromagnetism.jl")

#Plotting and graphics
include("plot/plot.jl")

# Precompile hints
include("precompile.jl")

end # module
