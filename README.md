# MultipleScattering

[![DOI](https://zenodo.org/badge/96763392.svg)](https://zenodo.org/badge/latestdoi/96763392)
[![Build Status](https://travis-ci.org/jondea/MultipleScattering.jl.svg?branch=master)](https://travis-ci.org/jondea/MultipleScattering.jl)
[![Coverage Status](https://coveralls.io/repos/github/jondea/MultipleScattering.jl/badge.svg?branch=master)](https://coveralls.io/github/jondea/MultipleScattering.jl?branch=master)
[![codecov.io](http://codecov.io/github/jondea/MultipleScattering.jl/coverage.svg?branch=master)](http://codecov.io/github/jondea/MultipleScattering.jl?branch=master)

A Julia library for simulating, processing and plotting multiple scattering of waves. 

The library uses the multipole method to solve the Helmholtz equation
(time-harmonic acoustics) in two dimensions.
In short, the method solves the problem with a series of Hankel functions
positioned at each particle, where the coefficents are picked so that the
boundary conditions are satisfied on the particle boundaries.
For a more lucid and complete explanation, see [Martin (1995)](https://pdfs.semanticscholar.org/8bd3/38ec62affc5c89592a9d6d13f1ee6a7d7e53.pdf).

This method is particularly efficient at solving acoustic problems with lots of
circular scatterers set in an infinite domain.

## Get started
This package is tested and works for Julia 0.6 and 0.5.
To get started, download and include the library
```julia
Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
using MultipleScattering
```

## Simple example
### Run
Define two particles, the first centred at [-2.,2.] and the second at [-2.,-2.]
```julia
p1 = Particle([-2.,2.])
p2 = Particle([-2.,-2.])
particles = [p1,p2]
```

Specify the angular frequency of the defualt incident plane wave ![incident plane wave](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20e%5E%7Bi%20%28k%20x%20-%20%5Comega%20t%29%7D) and calculate the response
```julia
w_arr = collect(0.1:0.01:1.)
simulation = FrequencySimulation(particles, w_arr; source_direction = [1.0,0.0])
```

### Plot
The package also provides recipes to be used with the `Plots` package for
plotting simulations after they have been run.
In our above simulation we ran the simulation for 100 different wavenumbers, and
measured the response at the default location (-10.,0.).
We can plot the time-harmonic response across these wavenumbers by typing:
```julia
using Plots
pyplot()
plot(simulation)
```
![Plot of response against wavenumber](example/intro/plot_simulation.png)

For a better overview you can plot the whole field for a specific k by typing:
```julia
plot(simulation,0.8)
```
![Plot real part of acoustic field](example/intro/plot_field.png)

This measures the field at lots of points in the domain rather than just at the listener position (the big green ball).
This way we can get an understanding of what is happening for one particular
wavenumber.

Note: most things in the package can be plotted by typing `plot(thing)` if you
need an insight into a specific part of your simulation.

## More examples
There are a lot of defaults implicit in this basic example.
Almost every part of the problem can be controlled, for example we can manually
construct the set of particles, define their positions, radii and give them
specific material properties. For all examples see [here](example/README.md).

## Acknowledgements and contributing
This library was restructured from one written by Artur Gower and Jonathan
Deakin.
Please contribute, if nothing else, criticism is welcome.
We are relatively new to Julia, and this is our first package, if anything is
untoward or even non-standard, please let us know.
