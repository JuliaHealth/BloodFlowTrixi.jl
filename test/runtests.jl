using BloodFlowTrixi
using Test

@testset "BloodFlowTrixi.jl" begin
    include("boundary_conditions.jl")
    include("model_equations.jl")
    include("./Aqua/aquatest.jl")
    include("./Extensions/DataInterpolationsTest.jl")

    @testset "1D Blood Flow Model" begin
        include("../exemples/Model1D/exemple.jl")
        @test sol.t[end] == tspan[end]
        @test all(u -> all(isfinite, u), sol.u)
    end
    @testset "1D order 2 Blood Flow Model" begin
        include("../exemples/Model1DOrd2/exemple.jl")
        @test sol.t[end] == tspan[end]
        @test all(u -> all(isfinite, u), sol.u)
    end
    @testset "2D Blood Flow Model" begin
        include("../exemples/Model2D/exemple.jl")
        @test sol.t[end] == tspan[end]
        @test all(u -> all(isfinite, u), sol.u)
    end
end
