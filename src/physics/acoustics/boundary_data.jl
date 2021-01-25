function boundary_data(shape1::Shape{T}, inside_medium::Acoustic{T,Dim}, outside_medium::Acoustic{T,Dim},
        sim::FrequencySimulation{T,Dim,P}, ωs::Vector{T};
        dr = T(10)^(3*Dim) * eps(T), kws...) where {T<:AbstractFloat,Dim,P<:Acoustic{T,Dim}}

    in_m = inside_medium
    out_m = outside_medium

    # points just inside particles
    inside1_points = boundary_points(shape1; dr = - dr - 10*eps(T))[:]
    inside2_points = boundary_points(shape1; dr = - 10*eps(T))[:]
    inside_points = mean([inside1_points,inside2_points])

    # points just outside particles
    outside1_points = boundary_points(shape1; dr = 10*eps(T))[:]
    outside2_points = boundary_points(shape1; dr = dr + 10*eps(T))[:]
    outside_points  = mean([outside1_points,outside2_points])

    # results from just inside particles
    L = length(inside1_points)
    N = 6; #number of sets of points
    all_points = [inside1_points; inside2_points; inside_points; outside1_points; outside2_points; outside_points]

    results = run(sim, all_points, ωs; kws...)
    in1 = field(results)[1:L,:]
    in2 = field(results)[L+1:2L,:]
    in  = field(results)[2L+1:3L,:]

    out1 = field(results)[3L+1:4L,:]
    out2 = field(results)[4L+1:5L,:]
    out  = field(results)[5L+1:6L,:]

    in_pressure = FrequencySimulationResult(in, inside_points, Vector(ωs))
    out_pressure = FrequencySimulationResult(out, outside_points, Vector(ωs))

    fields = (in2 - in1)/(dr * in_m.ρ)
    in_displace = FrequencySimulationResult(fields, (inside1_points + inside2_points) ./ 2, Vector(ωs))

    fields = (out2 - out1)/(dr * out_m.ρ)
    out_displace = FrequencySimulationResult(fields, (outside1_points + outside2_points) ./ 2, Vector(ωs))

    return ([in_pressure, out_pressure], [in_displace, out_displace])
end
