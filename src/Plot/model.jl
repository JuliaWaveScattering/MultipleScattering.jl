
"Plot the response across all wavenumbers"
function plot_model(model::FrequencyModel)
    plot(model.k_arr, [real(model.response) imag(model.response)], label=["real" "imaginary"], xlabel="Wavenumber (k)", ylabel="Response")
    plot!(title="a=$(model.particles[1].r), volfrac=$(calculate_volfrac(model)), shape=$(name(model.shape)) at ($(model.listener_positions[1,1]), $(model.listener_positions[2,1]))")
end