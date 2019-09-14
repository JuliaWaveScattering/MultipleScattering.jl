# Time response

```@meta
DocTestSetup = quote
    using MultipleScattering
end
```

This package calculates all scattering in the frequency domain, and we call the resulting field the frequency response $\hat u(\mathbf x,\omega)$, which satisfies $\nabla^2 \hat u(\mathbf x,\omega) + k^2 \hat u(\mathbf x,\omega) = 0$, where $k = \omega/c$. We can transform the frequency response $\hat u(\mathbf x,\omega)$ into a time response $u(\mathbf x,t)$ using a Fourier transform, where $u(\mathbf x,t)$ satisfies $\nabla^2 u(\mathbf x,t) - \frac{1}{c^2}  \frac{\partial^2}{\partial t^2}u(\mathbf x,t) = 0$. For a minimal example see [Results in time](@ref), or see [Technical details](@ref impulse_details) for more maths.

!!! note
    The package assumes the time response $u(\mathbf x,t)$ is always real, this simplifies the Fourier transform.

# Impulse function

As an example, let use a plane-wave source $\mathrm e^{\mathrm i \omega x}$ and measure the response at origin of the source $x = (0,0)$,
```jldoctest time
julia> plane_wave = plane_source(Acoustic(1.0, 1.0, 2); direction = [1.0,0.0], position = [0.0,0.0]);

julia> x = [[0.0,0.0]];

julia> ωs = LinRange(0.0,1.0,100);

julia> freq_response = run(plane_wave, x, ωs);

julia> t_vec = LinRange(0.0,100.,100);

julia> time_response = frequency_to_time(freq_response; t_vec = t_vec);
```
where we specified the times `t_vec` to calculate `time_response`. If no `t_vec` is given, the default times would be `t_vec = ω_to_t(ωs)` which is the standard choice for the Discrete Fourier Transform.  

Let us have a look at these responses:
```julia
julia> p1 = plot(freq_response, xlims = (0,2), ylims = (0.,1.5), field_apply = real);

julia> p2 = plot(time_response, field_apply = real);
```
impulse = TimeDiracImpulse(0.0)

disc_impulse = DiscreteImpulse(t_vec, impulse.in_time.(t_vec), ωs, impulse.in_freq.(ωs))

plot(freq_response, x = [[0,0]],
      xlims = (0,2), ylims = (0.,1.5),
      field_apply = real)

When we calculate the time response, range of frequencies `ωs`,

#  [Technical details](@id impulse_details)

We can calculate the time response $u(t)$, from the frequency response $\hat u( \omega)$ by approximating the Fourier transform:

$u(t) = \frac{1}{2\pi} \int_{-\infty}^\infty \hat u(\omega)\mathrm e^{-\mathrm i \omega t} d\omega, \quad \hat u(\omega) = \int_{-\infty}^\infty u(t)\mathrm e^{\mathrm i \omega t} dt,$

where the second equation is the inverse transform.
In practice
