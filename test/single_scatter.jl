function single_scatter_test()
    print("Testing single scatterer against an analytic solution with a random position and two random radii across a range of wavenumbers...")

    listener = [-10.0, 0.0];
    xp = 10*rand(2);
    a1 = 10*rand();
    a2 = 10*rand();
    k_arr = linspace(0.01,2.0,200)
    M1 = Int(round(a1)); #largest hankelh1 order
    M2 = Int(round(a2)); #largest hankelh1 order

    # incident plane wave
    uin(x,k) = exp(im * k * (x[1] - listener[1]))

    # known solution for the scattering coefficients for a single scattering
    function exact_single_scatter(x, k, a, M)
        A(n,k) = -exp(im * k * (xp[1] - listener[1]) + im * n * pi / 2);
        θ = atan2(x[2] - xp[2], x[1] - xp[1])
        r = norm(x - xp)
        u = 0.0im
        for m = -M:M
            u += A(m, k) * besselj(m, a * k) / hankelh1(m, a * k) * hankelh1(m, k * r) * exp(im * m * θ)
        end
        return u
    end
    u1(x,k) = exact_single_scatter(x, k, a1, M1)
    u2(x,k) = exact_single_scatter(x, k, a2, M2)

    resp1 = [u1(listener, k) for k in k_arr];
    resp2 = [u2(listener, k) for k in k_arr];

    # Add on incoming wave for the total solution
    resp1 += [uin(listener, k) for k in k_arr];
    resp2 += [uin(listener, k) for k in k_arr];

    # Using the MultipleScattering package
    model1 = FrequencyModel([Particle(xp,a1)], collect(k_arr);
                            listener_positions=listener,
                            source_direction=[1.0, 0.0],
                            hankel_order=M1
    );

    # Build a similar model with a different radius and hankel_order
    model2 = FrequencyModel([Particle(xp,a2)], collect(k_arr);
                            listener_positions=listener,
                            source_direction=[1.;0.],
                            hankel_order=M2
    );

    # Take the L^2 error between the analytic solution and our model
    err1 = mean(abs.(resp1 - model1.response[:]));
    err2 = mean(abs.(resp2 - model2.response[:]));

    print("completed for radius=$a1 and $a2 with errors $err1 and $err2 \n")
    if err1 < 1e-9 && err2 < 1e-9
        return true
    else
        return false
    end

end
