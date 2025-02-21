using BloodFlowTrixi
using Test
using Aqua

@testset "BloodFlowTrixi.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(BloodFlowTrixi; ambiguities = false,)
    end
    @testset "1D Blood Flow Model" begin
        include("../exemples/Model1D/exemple.jl")
    end
    @testset "2D Blood Flow Model" begin
        include("../exemples/Model2D/exemple.jl")
    end
end
