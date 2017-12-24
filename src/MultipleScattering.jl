module MultipleScattering

export Shape, Circle, Rectangle, TimeOfFlight, volume, inside, bounding_box,
       name, boundary_functions,
       Particle, volume, array_of_particles, random_particles, random_particles!, isequal,
       mean_radius, std_radius,
       FrequencySimulation, calculate_volfrac, generate_responses!,
       StatisticalMoments,
       TimeSimulation, frequency_to_time, time_to_frequency, ω_to_t, t_to_ω, delta_freq_impulse, get_gaussian_freq_impulse


abstract Shape

import Base.isequal, Base.(==)
import SpecialFunctions: besselj, hankelh1

using RecipesBase
# using ProgressMeter

# Functions and types for defining the domains we solve on
include("domain/particle.jl")
include("domain/shape.jl")
include("domain/random_particles.jl")

# Type definition to hold our model
include("frequency_model.jl")
# Where the actual maths (multipole method) happens
include("response.jl")
# Different constructor definitions which run the whole thing
include("frequency_model_constructors.jl")

# Turn the frequency model into a time model
include("time_model.jl")
# Statistical analysis of a batch of models
include("moments.jl")

# Recipes for plotting
include("plot.jl")

end # module
