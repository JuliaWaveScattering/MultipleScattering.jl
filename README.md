# MultipleScattering

[![Build Status](https://travis-ci.org/jondea/MultipleScattering.jl.svg?branch=master)](https://travis-ci.org/jondea/MultipleScattering.jl)
[![Coverage Status](https://coveralls.io/repos/github/jondea/MultipleScattering.jl/badge.svg?branch=master)](https://coveralls.io/github/jondea/MultipleScattering.jl?branch=master)
[![codecov.io](http://codecov.io/github/jondea/MultipleScattering.jl/coverage.svg?branch=master)](http://codecov.io/github/jondea/MultipleScattering.jl?branch=master)

A Julia library for simulating, processing and plotting acoustic data from
scattering problems, particularly random ones.

The library uses the multipole method to solve the Helmholtz equation 
(time-harmonic acoustics) in two dimensions.
In short, the method solves the problem with a series of Hankel functions 
positioned at each particle, where the coefficents are picked so that the 
boundary conditions are satisfied on the particle boundaries.
For a more lucid and complete explanation, see [Martin (1995)](https://pdfs.semanticscholar.org/8bd3/38ec62affc5c89592a9d6d13f1ee6a7d7e53.pdf).

This method is particularly efficient at solving acoustic problems with lots of
circular scatterers set in an infinite domain.

## Get started
To get started, download and include the basic library
```julia
Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
using MultipleScattering
```
## Very basic example
### Run
The basic way to construct and compute a model is to specify the volume
fraction, particle radius and wavenumbers:
```julia
volfrac = 0.01
radius = 1.0
k_arr = collect(linspace(0.01,1.0,100))
model = FrequencyModel(volfrac,radius,k_arr)
```
### Plot
The package also provides recipes to be used with the `Plots` package for 
plotting models after they have been run.
In our above model we ran the simulation for 100 different wavenumbers, and 
measured the response at a location away from the side wall.
We can plot the time-harmonic response across these wavenumbers by typing:
```julia
using Plots
plot(model)
```

This is a rather abstract plot, for some context, you can plot the whole field
at a specific k by typing.
This runs an identical model, but measures the field at lots of points in the 
domain rather than just at the listener position.
This way we can get an understanding of what is happening for one particular 
wavenumber.

```julia
plot(model,0.3)
```

## Less basic example
However, there are a great deal of defaults implicit in this statement.
For finer grain control of the problem, we can construct our set of particles,
define their positions, radii and give them material properties
```julia
particles = Vector{Particle{Float64}}(2)
# Phase speed inside the particle
c = 1.0 + 0.0*im
# Density of particle
ρ = 1.0
# Define their positions and radii
particles[1] = Particle{Float64}([0.2,0.8],0.1,c,ρ)
particles[2] = Particle{Float64}([0.2,-0.8],0.2,c,ρ)
```

We can then run this model, and plot the resultant field:
```julia
two_particle_model = FrequencyModel(particles,k_arr)

plot(two_particle_model,0.5)
```

This library was restructured from one written by Artur Gower and Jonathan 
Deakin.
Please contribute, if nothing else, criticism is welcome.
We are relatively new to Julia, and this is our first package, if anything is
untoward or even non-standard, please let us know.
