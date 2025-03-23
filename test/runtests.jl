using BloodFlowTrixi
using Test

@testset "BloodFlowTrixi.jl"  begin
    include("./Aqua/aquatest.jl")
    include("./Extensions/DataInterpolationsTest.jl")

    @testset "1D Blood Flow Model" begin
        include("../exemples/Model1D/exemple.jl")
    end
    @testset "2D Blood Flow Model" begin
        include("../exemples/Model2D/exemple.jl")
    end
end