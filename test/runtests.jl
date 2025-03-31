using ProjectionOntoSLn
using Test
using Aqua
using ExplicitImports
using LinearAlgebra: norm

@testset "HyperbolicTransformation" begin

    ### test the matrix B ###
    B = ProjectionOntoSLn.B(3)
    @test B == [2 1; -1 1; -1 -2]

    ### test type limits of B ###
    try
        # this is too big for the Int8 index type
        ProjectionOntoSLn.B(130)

        # we should not get here
        @test false
    catch
        @test true
    end

    try
        # this is should work
        ProjectionOntoSLn.B(130, Int16)
        @test true
    catch
        @test false
    end
end

@testset "Projection methods" begin

    for method in [unconstrained_newton, constrained_newton, composite_step, root_finding]

        a = [2.0, 1.0, 0.5] # p = a !!
        result = project_onto_sln(method, a; maxIter = 100, tolerance = 1.0e-10, debug = false)

        @test a ≈ result.projection

        a = [4.335260478477755, 0.8314927421281216]
        result = project_onto_sln(method, a; maxIter = 100, tolerance = 1.0e-10, debug = false)
        p = result.projection

        @test prod(p) ≈ 1
        @test testStationarity(a, p)
    end
end

@testset "default methods" begin

    a = [2.0, 1.0, 0.5] # p = a !!
    result = project_onto_sln(a)

    @test a ≈ result.projection

    a = [4.335260478477755, 0.8314927421281216]
    result = project_onto_sln(a)
    p = result.projection


    @test prod(p) ≈ 1
    @test testStationarity(a, p)
end


@testset "initial value" begin

    a = [2.0, 1.0, 0.5] # p = a !!
    p = ProjectionOntoSLn.find_inital_value(a)
    @test a ≈ p

    for a in ([4, 0.8], [1.5, 1.5, 0], [100, 99, 0.1], [0.1, 0.01, 0.0001])
        p⁰ = ProjectionOntoSLn.find_inital_value(a)
        result = project_onto_sln(a)
        p = result.projection

        @test prod(p⁰) ≈ 1
        @test norm(p⁰ - p) / norm(p) < 1.0e-1 #should be close
    end
end


@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(ProjectionOntoSLn) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(ProjectionOntoSLn) === nothing
end

@testset "Aqua Quality Checks" begin
    Aqua.test_all(ProjectionOntoSLn)
end
