@testset "Shape" begin
    @testset "Rectangle" begin
        rectangle = Rectangle((0.0,0.0), (2.0,3.0))
        @test volume(rectangle) ≈ 6.0
        @test name(rectangle) == "Rectangle"

        @testset "Boundary functions" begin
            x, y = boundary_functions(rectangle)
            @test x.(0:0.1:1) ≈ [0.0, 0.8, 1.6, 2.0, 2.0, 2.0, 1.2, 0.4, 0.0, 0.0, 0.0]
            @test y.(0:0.1:1) ≈ [0.0, 0.0, 0.0, 0.6, 1.8, 3.0, 3.0, 3.0, 2.4, 1.2, 0.0]
            @test_throws(DomainError,x(-eps(Float64)))
            @test_throws(DomainError,x(1.0 + eps(Float64)))
            @test_throws(DomainError,y(-eps(Float64)))
            @test_throws(DomainError,y(1.0 + eps(Float64)))
        end
    end

    @testset "Circle" begin
        circle = Circle([6.7,8.9],2.0)
        circle_bounding_rectangle = bounding_rectangle(circle)
        @test volume(circle)/volume(circle_bounding_rectangle) ≈ 0.7853981633974483
        @test name(circle) == "Circle"

        @testset "Boundary functions" begin
            x, y = boundary_functions(Circle([-1.0,2.0],3.0))
            @test x.(0:0.1:1) ≈
            [2.0, 1.4270509831248424, -0.07294901687515765, -1.927050983124842, -3.427050983124842, -4.0, -3.4270509831248424, -1.9270509831248428, -0.07294901687515831, 1.427050983124842, 2.0]
            @test y.(0:0.1:1) ≈
            [2.0, 3.763355756877419, 4.853169548885461, 4.853169548885461, 3.7633557568774196, 2.0000000000000004, 0.2366442431225808, -0.8531695488854605, -0.8531695488854609, 0.23664424312257992, 1.9999999999999993]
            @test_throws(DomainError,x(-eps(Float64)))
            @test_throws(DomainError,x(1.0 + eps(Float64)))
            @test_throws(DomainError,y(-eps(Float64)))
            @test_throws(DomainError,y(1.0 + eps(Float64)))
        end
    end

    @testset "Time of flight" begin
        time_of_flight = TimeOfFlight([-10.0,0.0],40.0)
        time_of_flight_bounding_rectangle = bounding_rectangle(time_of_flight)
        ratio = volume(time_of_flight)/volume(time_of_flight_bounding_rectangle)
        # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
        @test ratio > 0.5
        @test ratio < 1.0
        @test name(time_of_flight) == "Time of flight from planar source"

        @testset "Boundary functions" begin
            # This a nice shape because the dimensions form a Pythagorean triple
            x, y = boundary_functions(TimeOfFlight([-3.0,0.0],8.0))
            @test x.(0:0.125:1) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.75, 0.0]
            @test y.(0:0.125:1) ≈ [-4.0, -2.0, 0.0, 2.0, 4.0, 2.0, 0.0, -2.0, -4.0]
            @test_throws(DomainError,x(-eps(Float64)))
            @test_throws(DomainError,x(1.0 + eps(Float64)))
            @test_throws(DomainError,y(-eps(Float64)))
            @test_throws(DomainError,y(1.0 + eps(Float64)))
        end
    end

    @testset "Time of flight from point" begin
        time_of_flight = TimeOfFlightFromPoint([-10.0,0.0],40.0)
        time_of_flight_bounding_rectangle = bounding_rectangle(time_of_flight)
        ratio = volume(time_of_flight)/volume(time_of_flight_bounding_rectangle)
        # Geometric arguments dictate that the ratio must be between 0.5 and 1.0
        @test ratio > 0.5
        @test ratio < 1.0
        @test name(time_of_flight) == "Time of flight from point source"

        @testset "Boundary functions" begin
            x, y = boundary_functions(TimeOfFlightFromPoint([-1.0,0.0],3.0))
            @test x.(0:0.1:1) ≈ [0.0,1.2182846930812632,1.9095426112571139,1.9095426112571139,1.2182846930812632,0.0,0.0,0.0,0.0,0.0,0.0]
            @test y.(0:0.1:1) ≈ [-2.82842712474619,-2.019706171808505,-0.7311373286046446,0.7311373286046446,2.0197061718085054,2.82842712474619,1.6970562748477145,0.5656854249492386,-0.5656854249492386,-1.6970562748477145,-2.8284271247461903]
            @test_throws(DomainError,x(-eps(Float64)))
            @test_throws(DomainError,x(1.0 + eps(Float64)))
            @test_throws(DomainError,y(-eps(Float64)))
            @test_throws(DomainError,y(1.0 + eps(Float64)))
        end
    end

    @testset "Plot Shapes" begin
        # Just try each to see if we have any errors (yes thats a very low bar)
        rectangle = Rectangle([0.0,0.0],[2.0,3.0])
        # plot(rectangle)

        circle = Circle([-1.0,2.0],2.0)
        # plot!(circle)

        timeofflight = TimeOfFlightFromPoint([-1.0,0.0],3.0)
        # plot!(timeofflight)

        @test true
    end
end
