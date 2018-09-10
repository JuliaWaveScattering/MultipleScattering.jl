# MultipleScattering

[![DOI](https://zenodo.org/badge/96763392.svg)](https://zenodo.org/badge/latestdoi/96763392)
[![Build Status](https://travis-ci.org/jondea/MultipleScattering.jl.svg?branch=master)](https://travis-ci.org/jondea/MultipleScattering.jl)
[![Coverage Status](https://coveralls.io/repos/github/jondea/MultipleScattering.jl/badge.svg?branch=master)](https://coveralls.io/github/jondea/MultipleScattering.jl?branch=master)
[![codecov.io](http://codecov.io/github/jondea/MultipleScattering.jl/coverage.svg?branch=master)](http://codecov.io/github/jondea/MultipleScattering.jl?branch=master)

A Julia library for simulating, processing, and plotting multiple scattering of acoustic waves.

The library uses the multipole method to solve the Helmholtz equation
(time-harmonic acoustics) in two dimensions. This method is particularly efficient at solving scattering problems for particles in an infinite domain. Currently the library is configured for circular particles with any radius, density, sound speed and packing fraction. For details on the maths see [Martin (1995)](https://pdfs.semanticscholar.org/8bd3/38ec62affc5c89592a9d6d13f1ee6a7d7e53.pdf) and [Gower et al. (2017)](https://arxiv.org/abs/1712.05427).

#### Near Surface Backscattering
If you are here to learn about
[Near Surface Backscattering](example/near_surface_backscattering), then [click here](example/near_surface_backscattering) to see an example. For details on the maths see [Gower et al. (2018)](https://arxiv.org/abs/1801.05490). To see how to take the [moments](example/moments) of the backscattering [click here](example/moments).

## Get started
This package is tested and works for Julia 0.6.
To get started, download and include the library
```julia
Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
using MultipleScattering
```

## Simple example
### Run
Define the properties of your host medium, for example
```julia
host_medium = Acoustic(1.0, 1.0, 2)
```
an acoustic medium in 2D with density 1 and wavespeed 1.

Next, define two dense, circular acoustic particles, the first centred at [-2,2] with radius 2 and the second at [-2,-2] with radius 0.5,
```julia
particle_medium = Acoustic(10.0, 1.0, 2)
p1 = Particle(particle_medium, Circle([-2.0,2.0], 2.0))
p2 = Particle(particle_medium, Circle([-2.0,-2.0], 0.5))
particles = [p1,p2]
```

Lastly we define the source, for example an incident plane wave (![incident plane wave](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20e%5E%7Bi%20%28k%20x%20-%20%5Comega%20t%29%7D)) using a helper function.
```julia
source = plane_source(host_medium, [0.0,0.0])
```

Once we have these three components, we can build our `FrequencySimulation` object
```julia
simulation = FrequencySimulation(host_medium, particles, source)
```

To get numerical results, we run our simulation for specific positions and angular frequencies,
```julia
x = [[-10.0,0.0], [0.0,0.0]]
ω = 0.01:0.01:1.0
result = run(simulation, x, ω)
```

### Plot
The package also provides recipes to be used with the `Plots` package for
plotting simulations after they have been run.
In our above simulation we ran the simulation for 100 different wavenumbers, and
measured the response at the location (-10,0).
We can plot the time-harmonic response across these wavenumbers by typing:
```julia
using Plots
plot(result)
```
![Plot of response against wavenumber](example/intro/plot_result.png)

For a better overview you can plot the whole field in space for a specific angular frequency by typing:
```julia
ω = 0.8
plot(simulation,ω)
```
![Plot real part of acoustic field](example/intro/plot_field.png)

This measures the field at lots of points in the domain, so we can get an
understanding of what is happening for one particular angular frequency.

Note: most things in the package can be plotted by typing `plot(thing)` if you
need an insight into a specific part of your simulation.

To calculate an incident plane wave pulse in time use:
```julia
time_result = frequency_to_time(result)
plot(time_result)
```
![Plot real part of acoustic field](example/intro/plot_time_result.png)

## More examples
There are a lot of defaults implicit in this basic example.
Almost every part of the problem can be controlled, for example we can manually
construct the set of particles, define their positions, radii and give them
specific material properties. For all examples see [here](example/README.md).

## Acknowledgements and contributing
This library was restructured from one written by [Artur L Gower](https://arturgower.github.io/) and Jonathan
Deakin.
Please contribute, if nothing else, criticism is welcome.
We are relatively new to Julia, and this is our first package, if anything is
untoward or even non-standard, please let us know.
