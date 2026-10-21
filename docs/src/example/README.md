# Overview

Each example below has its own folder in [docs/src/example](https://github.com/JuliaWaveScattering/MultipleScattering.jl/tree/master/docs/src/example) with the full code. Most examples need the package `Plots` for the figures.

## [Helmholtz resonator](helmholtz-resonator/resonator.md)
- Build a Helmholtz resonator from small particles placed in the shape of a bottle.
- Plot the field for different frequencies to find the resonant frequencies, where the field inside the cavity is large.

## [Simple random particles](random_particles/README.md)
- Place random particles in a rectangle, calculate the response at a listener for many frequencies, and plot the field.

## [Random particles in a circle](particles_in_circle/README.md)
- Compare the wave scattered by one big circle with the wave scattered by a circle filled with small particles.

## [Box size](box_size/README.md)
- Find which particles need to be simulated to calculate the response at a listener for a given time interval.

## [Convergence of the basis order](hankel_convergence/README.md)
- Check how fast the scattered wave converges when increasing the number of Hankel functions, that is, the basis order.

## [Statistical moments](moments/README.md)
- Extract statistical information from many simulations with random particles: the mean, standard deviation, skew and kurtosis (also known as moments).

## [Near surface backscattering](near_surface_backscattering/README.md)
- A method to calculate the backscattering from an infinite halfspace filled with particles.
- See backscattering in frequency and time.

## [Time response from a single particle](time_response_single_particle/README.md)
- Calculate the response in time from a plane wave pulse which is reflected by one particle.

## [Two particles](two_particles/README.md)
- Calculate and plot the field scattered by two particles with different properties.
