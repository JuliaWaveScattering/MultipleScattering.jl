__precompile__()

module MultipleScattering

export Shape, Circle, Rectangle, volume, name, congruent, inside,
       boundary_functions, boundary_points
export PhysicalProperties, Acoustic, Electromagnetic, AcousticCapsule,
       get_basis_function, dim, field_dim
export Particle, (==)
export Source, besselj_field, self_test, TwoDimAcousticPointSource, TwoDimAcousticPlanarSource,
       (*), (+)
export Simulation, TwoDimAcousticFrequencySimulation, run, TimeSimulation, forcing
export SimulationResult, FrequencySimulationResult, TimeSimulationResult
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
include("simulation.jl")
include("result.jl")
include("random/random.jl")
include("t_matrix.jl")
include("scattering_matrix.jl")

# Specific physical models
include("physics/diffbessel.jl")
include("physics/acoustics.jl")
include("physics/electromagnetism.jl")

# Precompile hints
include("precompile.jl")

end # module
