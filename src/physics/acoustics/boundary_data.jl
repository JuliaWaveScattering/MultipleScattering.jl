function boundary_data(shape::Shape{T,2}, inside_medium::Acoustic{T,2}, outside_medium::Acoustic{T,2},
        sim::FrequencySimulation{T,2,P}, ωs::Vector{T};
        dr = 1e6 * eps(T), kws...) where {T<:AbstractFloat, P<:Acoustic{T,2}}

    in_m = inside_medium
    out_m = outside_medium

    # points just inside particles
    inside1_points = boundary_points(shape; dr = - dr - 10*eps(T))
    inside2_points = boundary_points(shape; dr = - 10*eps(T))
    inside_points = mean([inside1_points,inside2_points])

    # points just outside particles
    outside1_points = boundary_points(shape; dr = 10*eps(T))
    outside2_points = boundary_points(shape; dr = dr + 10*eps(T))
    outside_points  = mean([outside1_points,outside2_points])

    # results from just inside particles
    in1_results = run(sim, inside1_points, ωs; kws...)
    in2_results = run(sim, inside2_points, ωs; kws...)
    in_pressure = run(sim, inside_points, ωs; kws...)

    fields = (in2_results.field - in1_results.field)/(dr * in_m.ρ)
    in_displace = FrequencySimulationResult(fields, inside_points, RowVector(ωs))

    # results from just outside particles
    out1_results = run(sim, outside1_points, ωs; kws...)
    out2_results = run(sim, outside2_points, ωs; kws...)
    out_pressure  = run(sim, outside_points, ωs; kws...)

    fields = (out2_results.field - out1_results.field)/(dr * out_m.ρ)
    out_displace = FrequencySimulationResult(fields, outside_points, RowVector(ωs))

    return ([in_pressure, out_pressure], [in_displace, out_displace])
end
