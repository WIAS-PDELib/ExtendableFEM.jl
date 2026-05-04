using Test
using ExtendableGrids
using ExtendableFEMBase
using ExtendableFEM
using ExampleJuggler
using ExplicitImports
using Metis
using Aqua
using Triangulate
using SimplexGridFactory
using Symbolics


include("test_dgblf.jl")
include("test_boundary_operator.jl")
include("test_itemintegrator.jl")
include("test_timedependence.jl")
include("test_nonlinear_operator.jl")
include("test_helper_functions.jl")

function run_examples()
    ExampleJuggler.verbose!(true)

    example_dir = joinpath(@__DIR__, "..", "examples")

    exclude_examples = [
        "Example204_LaplaceEVProblem.jl",
        "Example206_CoupledSubGridProblems.jl",
        "Example227_ObstacleProblemLVPP.jl",
        "Example264_StokesDarcy.jl",
        "Example282_IncompressibleMHD.jl",
        "Example284_LevelSetMethod.jl",
        "Example285_CahnHilliard.jl",
        "Example295_SlidingDroplet.jl",
        "Example330_HyperElasticity.jl",
    ]

    modules = filter(∉(exclude_examples), readdir(example_dir))
    return @testset "module examples" begin
        @testmodules(example_dir, modules)
    end
end

function run_all_tests()
    @testset "Aqua.jl" begin
        Aqua.test_all(
            ExtendableFEM;
            ambiguities = false,
            piracies = false,
        )
        Aqua.test_ambiguities(ExtendableFEM)
    end

    @testset "ExplicitImports" begin
        @test ExplicitImports.check_no_implicit_imports(ExtendableFEM) === nothing
        @test ExplicitImports.check_no_stale_explicit_imports(ExtendableFEM) === nothing
    end

    if isdefined(Docs, :undocumented_names) # >=1.11
        @testset "UndocumentedNames" begin
            @test isempty(Docs.undocumented_names(ExtendableFEM))
        end
    end


    run_boundary_operator_tests()
    run_dgblf_tests()
    run_nonlinear_operator_tests()
    run_itemintegrator_tests()
    run_dt_tests()
    run_test_helper_functions()

    return nothing
end

# run_all_tests()
run_examples()
