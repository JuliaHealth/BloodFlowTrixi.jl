using BloodFlowTrixi
using Test
using Aqua

@testset "BloodFlowTrixi.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(BloodFlowTrixi; ambiguities = false,)
    end
    # Write your tests here.
end
