module MultipleScattering

export Shape, Circle, Rectangle, TimeOfFlight, TimeOfFlightFromPoint, volume,
       inside, bounding_box, name, boundary_functions,
       Particle, volume, array_of_particles, random_particles, random_particles!, isequal,
       mean_radius, std_radius,
       FrequencySimulation, calculate_volfrac, generate_responses!,
       build_field_simulation,
       StatisticalMoments,
       TimeSimulation, frequency_to_time, time_to_frequency, ω_to_t, t_to_ω,
       delta_freq_impulse, get_gaussian_freq_impulse


abstract type Shape end

import Base.isequal, Base.(==)
import SpecialFunctions: besselj, hankelh1

using RecipesBase
using ProgressMeter

# Functions and types for defining the domains we solve on
include("domain/particle.jl")
include("domain/shape.jl")
include("domain/random_particles.jl")

# Type definition to hold our simulation
include("frequency_simulation.jl")
# Where the actual maths (multipole method) happens
include("response.jl")
# Different constructor definitions which run the whole thing
include("frequency_simulation_constructors.jl")

# Turn the frequency simulation into a time simulation
include("time_simulation.jl")
# Statistical analysis of a batch of simulations
include("moments.jl")

# Recipes for plotting
include("plot.jl")

end # module
