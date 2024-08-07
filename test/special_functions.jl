using GSL, SpecialFunctions
using Test, LinearAlgebra

@testset "Special functions" begin

    @testset "bessel functions" begin

        x = rand() - 0.5 + (rand() - 0.5) * im

        @test diffbessel(besselj,2,x,2) ≈ diffbesselj(2,x,2)
        @test diffbessel(hankelh1,2,x,2) ≈ diffhankelh1(2,x,2)

        @test diffbessel(besselj,4,x,3) ≈ diffbesselj(4,x,3)
        @test diffbessel(hankelh1,4,x,3) ≈ diffhankelh1(4,x,3)

        # spherical bessel functions
        @test shankelh1(1, x) ≈ - exp(im*x) * (x + im) / (x^2)
        @test sbesselj(1, x) ≈ sin(x)/x^2 - cos(x)/x

        @test sbesselj(1, eps(zero(Float64))) ≈ 0.0
        @test sbesselj(0, eps(zero(Float64))) ≈ 1.0

        n = 1
        @test 2 * diffsbessel(shankelh1,n,x) ≈ shankelh1(n-1, x) - (shankelh1(n, x) + x*shankelh1(n+1, x))/x

        @test diffsbessel(shankelh1,n,x) ≈ diffshankelh1(n,x)
        @test diffsbessel(sbesselj,n,x) ≈ diffsbesselj(n,x)

        n = 2
        @test 2 * diffsbessel(shankelh1,n,x) ≈ shankelh1(n-1, x) - (shankelh1(n, x) + x*shankelh1(n+1, x))/x

        xs = 0.01:0.9:11.0;
        zs = xs + (collect(xs) .* 1im)
        n = 15; ns = 0:n

        ε = eps(Float64) * 10^7

        maxerror = [
            abs(diffsbesselj(n,x) - (sbesselj(n,x+ε) - sbesselj(n,x-ε)) / (2ε)) / abs(diffsbesselj(n,x))
        for n in ns, x in zs] |> maximum    

        @test maxerror < 1e-4

        maxerror = [
            abs(diffshankelh1(n,x) - (shankelh1(n,x+ε) - shankelh1(n,x-ε)) / (2ε)) / abs(diffshankelh1(n,x))
        for n in ns, x in zs] |> maximum    

        @test maxerror < 1e-5

        maxerror = [
            abs(diff2sbesselj(n,x) - (diffsbesselj(n,x+ε) - diffsbesselj(n,x-ε)) / (2ε)) / abs(diff2sbesselj(n,x))
        for n in ns, x in zs] |> maximum    

        @test maxerror < 1e-4

        maxerror = [
            abs(diff2shankelh1(n,x) - (diffshankelh1(n,x+ε) - diffshankelh1(n,x-ε)) / (2ε)) / abs(diff2shankelh1(n,x))  
        for n in ns, x in zs] |> maximum    
    
        @test maxerror < 1e-5
        
        maxerror = [
            abs(diff3sbesselj(n,x) - (diff2sbesselj(n,x+ε) - diff2sbesselj(n,x-ε)) / (2ε)) / abs(diff3sbesselj(n,x))
        for n in ns, x in zs] |> maximum    

        @test maxerror < 1e-4

        maxerror = [
            abs(diff3shankelh1(n,x) - (diff2shankelh1(n,x+ε) - diff2shankelh1(n,x-ε)) / (2ε)) / abs(diff3shankelh1(n,x))  
        for n in ns, x in zs] |> maximum    
    
        @test maxerror < 1e-5

        # Compare with spherical bessel from GSL, which can only accept real arguments..
        errs = map(xs) do x
            maximum(abs.(sf_bessel_jl_array(n,x) - [sbesselj(n,x) for n in ns]))
        end;
        @test maximum(errs) < 20*eps(Float64)

        errs = map(xs) do x
            maximum(
                abs.(
                    im .* sf_bessel_yl_array(n,x) -
                    [shankelh1(n,x) for n in ns] +
                    sf_bessel_jl_array(n,x)
                ) ./ abs.(sf_bessel_yl_array(n,x))
            )
        end;
        @test maximum(errs) < 1e4*eps(Float64)

    end

    @testset "Associated legendre functions" begin
        x = rand(1)[1] * 0.99
        l_max = 2

        # a few associated legendre functions without the Condon-Shortly factor
        Plm_arr = [1,x,sqrt(1-x^2),(3x^2-1)/2, 3x*sqrt(1-x^2),3*(1-x^2)]

        ls, ms = associated_legendre_indices(l_max)

        @test sf_legendre_array(GSL_SF_LEGENDRE_NONE, l_max, x)[1:length(ls)] ≈ Plm_arr

        sph_factors = (-1).^ms .* sqrt.((2 .* ls .+ 1) ./ (4pi) .* factorial.(ls-ms) ./ factorial.(ls+ms))

        condon_phase = (-1).^ms

        @test condon_phase .* sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, x)[1:length(ls)] ≈ sph_factors .* Plm_arr

        # test for the complex case
        y = rand(1)[1] * 0.99
        z = x + y*1im

        # a few complex associated legendre functions
        Plm_arr_complex = [1,conj(sqrt(1-z^2))/2,z,-sqrt(1-z^2),conj(1-z^2)/8,conj(z*sqrt(1-z^2))/2,(3z^2-1)/2,-3z*sqrt(1-z^2),3*(1-z^2)]

        ls, ms = spherical_harmonics_indices(l_max)

        norm_factors = sqrt.((2 .* ls .+ 1) ./ (4pi) .* factorial.(ls - ms) ./ factorial.(ls + ms))

        @test complex_legendre_array(l_max, z) ≈ Plm_arr_complex .* norm_factors

    end

    @testset "Spherical harmonics" begin
        θ = rand(1)[1] * 0.99
        φ = rand(1)[1] * 0.99

        l_max = 8 # small l_max due to factorial formula below

        ls, ms = spherical_harmonics_indices(l_max)
        sphs = spherical_harmonics(l_max, θ, φ)
        sphs_dθ = spherical_harmonics_dθ(l_max, θ, φ)

        @test sphs_dθ[1] == Complex{Float64}(0)
        
        h = 1e5 * eps(Float64) 
        θ1 = θ - h
        sph1s = spherical_harmonics(l_max, θ1, φ)

        θ2 = θ + h
        sph2s = spherical_harmonics(l_max, θ2, φ)

        sphs_dθ_approx = (sph2s - sph1s) ./ (2h)

        @test maximum(abs.(sphs_dθ_approx - sphs_dθ)[2:end] ./ abs.(sphs_dθ[2:end])) < 5e-4
        
        @test maximum(i - lm_to_spherical_harmonic_index(ls[i],ms[i]) for i in eachindex(ls)) == 0

        # check special case l == abs(m)
        is = findall(ls .== abs.(ms))

        for i in is
            @test sphs[i] ≈ (- sign(ms[i]))^ls[i] / (2^ls[i] * factorial(ls[i])) *
                sqrt(factorial(2*ls[i] + 1) / (4pi)) * sin(θ)^ls[i] * exp(im * ms[i] * φ)
        end

        l_max = 30
        ls, ms = spherical_harmonics_indices(l_max)
        sphs = spherical_harmonics(l_max, θ, φ)

        # special case m == 0, reduce to just Legendre polynomials
        Ps = sf_legendre_Pl_array(l_max, cos(θ))
        is = findall(ms .== 0)

        for l in 0:l_max
            @test sphs[is][l+1] ≈ sqrt((2l+1)/(4pi)) * Ps[l+1]
        end

        #sphs[inds] .≈ sqrt.((2 .* (0:l_max) .+ 1) ./ (4pi)) .* Ps

        # special case, north pole
        θ = 0.0
        sphs = spherical_harmonics(l_max, θ, φ)

        for i in eachindex(sphs)
            if ms[i] == 0
                @test sphs[i] ≈ sqrt((2*ls[i] + 1)/(4pi))
            else
                @test sphs[i] ≈ 0.0
            end
            i += 1
        end

        # Test complex argument spherical harmonics
        sphs_complex = spherical_harmonics(l_max, θ+0.0im, φ)
        @test sphs ≈ sphs_complex

    end

    @testset "Gaunt and Wigner symbols" begin
        l1 = rand(1:100)
        l2 = rand(1:100)
        l3 = rand(abs(l1-l2):2:(l1+l2)) # guarantees that iseven(l1+l2+l3)

        @test iseven(l1+l2+l3)

        m1 = rand(-l1:l1)
        m2 = rand(-l2:l2)
        m3 = m1-m2

        @test_throws(DomainError,gaunt_coefficient(l1,2*l1,l2,m2,l3,m3))
        @test_throws(MethodError,gaunt_coefficient(l1,m1,l2,m2,0.1,m3))

        # the spherical harmonics linearisation formula
        θ, φ = rand(2) .* pi

        l_small = 6
        l_max = 2*l_small # needs to be larger than l_small

        ls, ms = spherical_harmonics_indices(l_max);
        Ys = spherical_harmonics(l_max, θ, φ);

        cs = reshape(
            [
                gaunt_coefficient(l1,m1,l2,m2,l3,m3)

                    for l3 = 0:l_max for m3 = -l3:l3
                for l2 = 0:l_small for m2 = -l2:l2
            for l1 = 0:l_small for m1 = -l1:l1]
        , ((l_small + 1)^2,(l_small + 1)^2,(l_max+1)^2));

        for n2 in 1:(l_small + 1)^2, n3 in 1:(l_small + 1)^2
            @test 4pi * (-1)^(ms[n3]+ms[n2]) * (1.0im)^(ls[n3]-ls[n2]) * Ys[n3] * conj(Ys[n2]) ≈
                sum( (1.0im).^(-ls) .* (-1.0).^ms .* conj.(Ys) .* cs[n2,n3,:]) atol = 1e-12
        end

    end

    @testset "complex radial coordinates" begin

    # Test 3-dimensional transforms
        xs = [rand(-1.01:0.1:1.0,3) + rand(-1.01:0.1:1.0,3)*im for i = 1:100]
        rθφs = cartesian_to_radial_coordinates.(xs)
        @test maximum(norm.(xs - radial_to_cartesian_coordinates.(rθφs))) < 2e-14
        
        rθφs = cartesian_to_spherical_coordinates.(xs)
        vs = [rand(-1.01:0.1:1.0,3) + rand(-1.01:0.1:1.0,3)*im for i = 1:100]

        svs = [
            cartesian_to_spherical_vector(vs[i],xs[i])
        for i in eachindex(vs)]
            
        v2s = [
            spherical_to_cartesian_vector(svs[i],rθφs[i])
        for i in eachindex(vs)]

        @test maximum(norm.(vs - v2s)) < 5e-14
        
        # In a Cartesian representation
        rθφ = cartesian_to_spherical_coordinates(xs[1])
        r, θ, φ = rθφ
        
        er = [cos(φ) * sin(θ), sin(φ) * sin(θ), cos(θ)]
        eθ = [cos(φ) * cos(θ), sin(φ) * cos(θ), -sin(θ)]
        eϕ = [-sin(φ), cos(φ), 0.0]

        @test [1.0,0.0,0.0] - cartesian_to_spherical_vector(er,xs[1]) |> norm < 1e-14
        @test [0.0,1.0,0.0] - cartesian_to_spherical_vector(eθ,xs[1]) |> norm < 1e-14
        @test [-.0,0.0,1.0] - cartesian_to_spherical_vector(eϕ,xs[1]) |> norm < 1e-14

        xs = [rand(-1.01:0.1:1.0,3) for i = 1:100]
        rθφs = cartesian_to_radial_coordinates.(xs)

        @test maximum(rθφ[2] for rθφ in rθφs) <= pi
        @test minimum(rθφ[2] for rθφ in rθφs) >= 0.0

        @test pi/2 < maximum(rθφ[3] for rθφ in rθφs) <= pi
        @test -pi <= minimum(rθφ[3] for rθφ in rθφs) < -pi/2


    # Test 2-dimensional transforms
        xs = [rand(-1.01:0.1:1.0,2) + rand(-1.01:0.1:1.0,2)*im for i = 1:100]
        rθs = cartesian_to_radial_coordinates.(xs)

        @test maximum(norm.(xs - radial_to_cartesian_coordinates.(rθs))) < 1e-14
    end
end
