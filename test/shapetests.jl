@testset "Shape" begin
    @testset "Rectangle" begin
        rectangle = Rectangle([0.0,0.0],[2.0,3.0])
        @test volume(rectangle) ≈ 6.0
        @test name(rectangle) == "Rectangle"

        @testset "Boundary functions" begin
            x, y = boundary_functions(rectangle)
            @test x.(0:0.1:1) ≈ [0.0, 0.8, 1.6, 2.0, 2.0, 2.0, 1.2, 0.4, 0.0, 0.0, 0.0]
            @test y.(0:0.1:1) ≈ [0.0, 0.0, 0.0, 0.6, 1.8, 3.0, 3.0, 3.0, 2.4, 1.2, 0.0]
            @test_throws(Exception,x(-eps(Float64)))
            @test_throws(Exception,x(1.0 + eps(Float64)))
            @test_throws(Exception,y(-eps(Float64)))
            @test_throws(Exception,y(1.0 + eps(Float64)))
        end
    end

    @testset "Circle" begin
        circle = Circle(2.0,[6.7,8.9])
        circle_bounding_box = bounding_box(circle)
        @test volume(circle)/volume(circle_bounding_box) ≈ 0.7853981633974483
        @test name(circle) == "Circle"

        @testset "Boundary functions" begin
            x, y = boundary_functions(Circle(3.0,[-1.0,2.0]))
            @test x.(0:0.1:1) ≈ [0.0,-0.19098300562505255,-0.6909830056250525,-1.3090169943749475,-1.8090169943749475,-2.0,-1.8090169943749475,-1.3090169943749475,-0.6909830056250528,-0.19098300562505266,0.0]
            @test y.(0:0.1:1) ≈ [2.0,2.5877852522924734,2.9510565162951536,2.9510565162951536,2.5877852522924734,2.0,1.412214747707527,1.0489434837048464,1.0489434837048464,1.4122147477075266,1.9999999999999998]
            @test_throws(Exception,x(-eps(Float64)))
            @test_throws(Exception,x(1.0 + eps(Float64)))
            @test_throws(Exception,y(-eps(Float64)))
            @test_throws(Exception,y(1.0 + eps(Float64)))
        end
    end

    @testset "Time of flight" begin
        time_of_flight = TimeOfFlight([-10.0,0.0],40.0)
        time_of_flight_bounding_box = bounding_box(time_of_flight)
        ratio = volume(time_of_flight)/volume(time_of_flight_bounding_box)
        # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
        @test ratio > 0.5
        @test ratio < 1.0
        @test name(time_of_flight) == "Time of flight"

        @testset "Boundary functions" begin
            # This a nice shape because the dimensions form a Pythagorean triple
            x, y = boundary_functions(TimeOfFlight([-3.0,0.0],8.0))
            @test x.(0:0.125:1) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.75, 0.0]
            @test y.(0:0.125:1) ≈ [-4.0, -2.0, 0.0, 2.0, 4.0, 2.0, 0.0, -2.0, -4.0]
            @test_throws(Exception,x(-eps(Float64)))
            @test_throws(Exception,x(1.0 + eps(Float64)))
            @test_throws(Exception,y(-eps(Float64)))
            @test_throws(Exception,y(1.0 + eps(Float64)))
        end
    end

    @testset "Time of flight from point" begin
        time_of_flight = TimeOfFlightFromPoint([-10.0,0.0],40.0)
        time_of_flight_bounding_box = bounding_box(time_of_flight)
        ratio = volume(time_of_flight)/volume(time_of_flight_bounding_box)
        # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
        @test ratio > 0.5
        @test ratio < 1.0
        @test name(time_of_flight) == "Time of flight from point"

        @testset "Boundary functions" begin
            x, y = boundary_functions(TimeOfFlightFromPoint([-1.0,0.0],3.0))
            @test x.(0:0.1:1) ≈ [0.0,1.2182846930812632,1.9095426112571139,1.9095426112571139,1.2182846930812632,0.0,0.0,0.0,0.0,0.0,0.0]
            @test y.(0:0.1:1) ≈ [-2.82842712474619,-2.019706171808505,-0.7311373286046446,0.7311373286046446,2.0197061718085054,2.82842712474619,1.6970562748477145,0.5656854249492386,-0.5656854249492386,-1.6970562748477145,-2.8284271247461903]
            @test_throws(Exception,x(-eps(Float64)))
            @test_throws(Exception,x(1.0 + eps(Float64)))
            @test_throws(Exception,y(-eps(Float64)))
            @test_throws(Exception,y(1.0 + eps(Float64)))
        end
    end

    @testset "Plot Shapes" begin
        # Just try each to see if we have any errors (yes thats a very low bar)
        rectangle = Rectangle([0.0,0.0],[2.0,3.0])
        plot(rectangle)

        circle = Circle(3.0,[-1.0,2.0])
        plot!(circle)

        timeofflight = TimeOfFlightFromPoint([-1.0,0.0],3.0)
        plot!(timeofflight)

        @test true
    end
end
