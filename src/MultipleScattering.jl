module MultipleScattering

export Shape, Circle, Rectangle, volume, name
export PhysicalProperties, Acoustic, Electromagnetic, AcousticCapsule,
       basis_function, dim, field_dim
export Particle, (==)
export Source, self_test, TwoDimAcousticPointSource, TwoDimAcousticPlanarSource,
       (*), (+)
export Simulation, TwoDimAcousticFrequencySimulation, run, TimeSimulation
export SimulationResult, FrequencySimulationResult, TimeSimulationResult
export SimulationDistribution, FrequencySimulationDistribution, sample


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
include("random.jl")

end # module
