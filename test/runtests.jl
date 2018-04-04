import Base.Test: @testset, @test, @test_throws

import StaticArrays: MVector

using MultipleScattering

@testset "Tests" begin
    x = MVector(1.0,1.0)
    c = Circle(2.0)

    # 2D Acoustics
    a2 = Acoustics(1.0,1.0 + 0.0im,2)
    @test dim(a2) == 2
    @test field_dim(a2) == 1

    # 3D Acoustics
    a3 = Acoustics(1.0,1.0 + 0.0im,3)
    @test dim(a3) == 3
    @test field_dim(a3) == 1

    p = Particle(x,a2,c)

    @test_throws MethodError Particle(x,a3,c)

    source_position = MVector(0.0,1.0)
    amplitude = 1.0
    s1 = TwoDimAcousticPointSource(a2, source_position, amplitude)
    s2 = TwoDimAcousticPointSource(a2, 2.*source_position, amplitude)

    s3 = 2*s1 + s2
end
