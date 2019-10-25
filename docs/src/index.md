# MultipleScattering.jl Documentation

*A Julia library for simulating, processing, and plotting multiple scattering of acoustic waves.*

The library uses the multipole method to solve the Helmholtz equation
(time-harmonic waves). The multipole method is particularly efficient at solving scattering problems for particles in an infinite domain. This library is configured to use T-matrices to represent scattering from particles with any shape and properties. The package is setup to deal with different spatial dimensions and types of waves which satisfy Helmholtz equation's, e.g. acoustics, electromagnetism, elasticity. For details on some of the maths see [Martin (1995)](https://pdfs.semanticscholar.org/8bd3/38ec62affc5c89592a9d6d13f1ee6a7d7e53.pdf) and [Gower et al. (2017)](https://arxiv.org/abs/1712.05427).

## Installation

Install Julia v1.0 or later, then run

```julia
using Pkg
Pkg.clone("https://github.com/jondea/MultipleScattering.jl.git")
using MultipleScattering
```

## Manual

You can learn to use this package through [examples](example/README.md) or through our manual, which starts with a simple [Introduction](@ref).


## Contents
```@contents
Pages = ["base.md", "acoustics.md","random.md"]
Depth = 2
```
