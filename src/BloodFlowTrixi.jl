# should you ask why the last line of the docstring looks like that:
# it will show the package path when help on the package is invoked like     help?> BloodFlowTrixi
# but will interpolate to an empty string on CI server, preventing appearing the path in the documentation built there

"""
    Package BloodFlowTrixi v$(pkgversion(BloodFlowTrixi))

This package implements 1D and 2D blood flow models for arterial circulation using Trixi.jl, enabling efficient numerical simulation and analysis.

Docs under https://yolhan83.github.io/BloodFlowTrixi.jl

$(isnothing(get(ENV, "CI", nothing)) ? ("\n" * "Package local path: " * pathof(BloodFlowTrixi)) : "") 
"""
module BloodFlowTrixi
    export flux_nonconservative,BloodFlowEquations1D,BloodFlowEquations2D,boundary_condition_outflow,boundary_condition_pressure_in,initial_condition_simple,source_term_simple
    using Trixi
    # Write your package code here.
    abstract type AbstractBloodFlowEquations{NDIMS, NVARS} <:Trixi.AbstractEquations{NDIMS, NVARS} end
    include("1DModel/1dmodel.jl")
    include("./1DModel/Test_Cases/pressure_in.jl")
    include("1DModel/Test_Cases/convergence_test.jl")
    include("2DModel/2dmodel.jl")
end
