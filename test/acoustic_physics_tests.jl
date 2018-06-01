
@testset "acoustic sources" begin
    a2 = Acoustic(0.1,0.1 + 0.0im,2)
    a2_host = Acoustic(1.0,1.0 + 0.0im,2)

    # Create two point sources
    source_position = SVector(0.0,1.0)
    amplitude = 1.0
    s1 = TwoDimAcousticPointSource(a2, source_position, amplitude)
    s2 = TwoDimAcousticPointSource(a2, 2.*source_position, amplitude)

    # Create new souce as a linear combination of two other sources
    s3 = 2*s1 + s2

    # Check that the field is indeed a linear conbination
    x = SVector(1.0, 1.0)
    @test s3.field(x,1.0) == 2*s1.field(x,1.0) + s2.field(x,1.0)

    # Test the bessel expansions of the source
    ω = 0.8
    centre =  SVector(1.0,0.0)
    s3_besselj = besselj_field(s3, a2, centre; basis_order = 7)
    xs = [centre + 0.1.*[cos(τ),sin(τ)] for τ = 0.0:0.3:1.5]
    @test norm([s3.field(x,ω) - s3_besselj(x,ω) for x in xs]) < 1e-7*norm([s3.field(x,ω) for x in xs])

    source = TwoDimAcousticPlanarSource(a2_host, SVector(-10.0,0.0), SVector(1.0,0.0), 1.0)
    source_besselj = besselj_field(source, a2_host, centre)
    @test norm([source.field(x,ω) - source_besselj(x,ω) for x in xs]) < 2e-9*norm([source.field(x,ω) for x in xs])
end
