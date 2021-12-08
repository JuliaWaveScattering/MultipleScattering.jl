import Test: @testset, @test, @test_throws
import Statistics: mean

using MultipleScattering
using LinearAlgebra

@testset "Fourier transform" begin

    ω_arr = collect(LinRange(0.0,100.0,200))
    t_arr = ω_to_t(ω_arr)
    N = length(ω_arr) - 1
    t0 = t_arr[N]

    # test invertablity
    @test ω_arr ≈ t_to_ω(ω_to_t(ω_arr))
    fs = rand(-1:0.01:1.,2N+1,1)
    @test fs ≈ frequency_to_time(time_to_frequency(fs,t_arr),ω_arr)

    # Gaussian to Gaussian
    as = [0.5];
    Fs = [exp(-a*ω^2)*exp(im*t0*ω) for ω in ω_arr, a in as]
    fs = frequency_to_time(Fs,ω_arr,t_arr);
    true_fs = (1/(2pi))*[sqrt(pi/a)*exp(-(t-t0)^2/(4*a)) for t in t_arr, a in as]
    @test fs ≈ true_fs
    @test Fs ≈ time_to_frequency(fs,t_arr,ω_arr)

    # sinc to rectangle
    rect(t) = (abs(t)<0.5) ? 1.0 : 0.
    as = [0.7]; #only positive values
    Fs = [ (ω == 0.0im) ? 1.0+0.0im : sin(a*pi*ω)*exp(im*t0*ω)/(a*pi*ω) for ω in ω_arr, a in as]
    fs = frequency_to_time(Fs,ω_arr,t_arr) # only correct for half of time_arr
    true_fs = (1/(2pi))*[(1/abs(a))*rect((t-t0)/(2*pi*a)) for t in t_arr, a in as]
    @test mean(abs.(fs-true_fs)) < 0.02*mean(abs.(true_fs))

    # trapezoidal integration
    ω_arr = LinRange(0.0,130.0,600)
    t_arr = ω_to_t(ω_arr)

    fs = cos.(t_arr).*exp.(-(t_arr .- t_arr[end]/2).^2/(t_arr[end]^2/25))
    fs = reshape(fs,length(fs),1)

    Fs_trap = time_to_frequency(fs, t_arr; method=:trapezoidal)
    Fs = time_to_frequency(fs, t_arr)
    # plot(ω_arr,[real(Fs),real(Fs_trap)])
    @test mean(abs.(Fs-Fs_trap)) < 0.0012*mean(abs.(Fs))

    fs_trap = frequency_to_time(Fs_trap,ω_arr; method=:trapezoidal)
    @test mean(abs.(fs-fs_trap))< 1e-5*mean(abs.(fs))

    ## rectangle to sinc
    # as = [1/(2*(maximum(ω_arr)+ω_arr[2]))]; #only positive values
    # # t0 = t_arr[end]/3
    # Fs = Complex{Float64}[ rect(a*ω)*exp(im*t0*ω) for ω in ω_arr, a in as]
    # fs =  frequency_to_time(Fs, ω_arr, t_arr)
    # true_fs = [ (t==t0) ? 1/(2pi*a) : sin((t-t0)/(2a))/(pi*(t-t0)) for t in t_arr, a in as]
    # @test maximum(abs.(fs-true_fs)) < 1e-3*maximum(abs.(true_fs))

    # pole to wavy heaviside
    # ω_arr = collect(LinRange(0.0,10.0,4000)) # need fine mesh near ω=0
    # dω = (ω_arr[2]-ω_arr[1])
    # N = length(ω_arr)-1;
    # time_arr = LinRange(0,2π/dω,2N+2)[1:(2N+1)]
    # as = [0.1]; #only positive values
    # Fs = Complex{Float64}[1/(a-im*ω) for ω in ω_arr, a in as]
    # fs =  frequency_to_time(Fs,ω_arr,time_arr);
    # true_fs = [exp(-a*t) for t in time_arr, a in as]
    # mean(abs.(fs-true_fs))
    #  1e-3*mean(abs.(true_fs))

end
