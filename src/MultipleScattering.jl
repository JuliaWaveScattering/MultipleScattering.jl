module MultipleScattering

export Shape, Circle, Rectangle, volume, name
export PhysicalProperties, Acoustic, Electromagnetic, AcousticCapsule,
       get_basis_function, dim, field_dim
export Particle, (==)
export Source, self_test, TwoDimAcousticPointSource, TwoDimAcousticPlanarSource,
       (*), (+)
export Simulation, TwoDimAcousticFrequencySimulation, run, TimeSimulation
export SimulationResult, FrequencySimulationResult, TimeSimulationResult
export SimulationDistribution, FrequencySimulationDistribution, sample
export t_matrix, get_t_matrices
export scattering_matrix

import SpecialFunctions: besselj, hankelh1

import StaticArrays: SVector

import OffsetArrays: OffsetArray

using RecipesBase
using ProgressMeter

include("shape.jl")
include("physics.jl")
include("particle.jl")
include("source.jl")
include("simulation.jl")
include("result.jl")
include("random.jl")
include("diffbessel.jl")
include("t_matrix.jl")
include("scattering_matrix.jl")

end # module
