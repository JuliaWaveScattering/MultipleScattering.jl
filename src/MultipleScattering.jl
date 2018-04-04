module MultipleScattering

export Shape, Circle, Rectangle, volume, name
export PhysicalProperties, Acoustics, Electromagnetism, AcousticCapsule,
       basis_function, dim, field_dim
export Particle, (==)
export Source, self_test, TwoDimAcousticPointSource, TwoDimAcousticPlanarSource,
       (*), (+)
export Simulation, TwoDimAcousticFrequencySimulation, run, TimeSimulation
export SimulationResult, FrequencySimulationResult, TimeSimulationResult

import SpecialFunctions: besselj, hankelh1

import StaticArrays: MVector

using RecipesBase
using ProgressMeter

include("shape.jl")
include("physics.jl")
include("particle.jl")
include("source.jl")
include("simulation.jl")
include("result.jl")

end # module
