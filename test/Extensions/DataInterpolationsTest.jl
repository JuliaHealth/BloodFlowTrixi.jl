using DataInterpolations

@testset "2D Blood Flow Model with interpolation" begin
    include("../../exemples/Model2D/diexemple.jl")
    @test sol.t[end] == tspan[end]
    @test all(u -> all(isfinite, u), sol.u)
end
