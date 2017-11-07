
"Plot the response across all wavenumbers"
function plot_model(model::FrequencyModel)
    plot(model.k_arr, [real(model.response) imag(model.response)], label=["real" "imaginary"], xlabel="Wavenumber (k)", ylabel="Response", grid=false)
    plot!(title="a=$(model.particles[1].r), volfrac=$(calculate_volfrac(model)), shape=$(name(model.shape)) at ($(model.listener_positions[1,1]), $(model.listener_positions[2,1]))")
end


"Plot the response across time"
function plot_model(model::TimeModel)
    plot(model.time_arr, [real(model.response) imag(model.response)], label=["real" "imaginary"], xlabel="Time (t)", ylabel="Response", grid=false)
    plot!(title="a=$(model.frequency_model.particles[1].r), volfrac=$(calculate_volfrac(model.frequency_model)), shape=$(name(model.frequency_model.shape)) at ($(model.frequency_model.listener_positions[1,1]), $(model.frequency_model.listener_positions[2,1]))")
end
