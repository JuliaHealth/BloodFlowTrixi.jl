using Trixi
using OrdinaryDiffEq
using DataInterpolations
using BloodFlowTrixi
using StaticArrays, LinearAlgebra
using QuadGK
using LinearAlgebra

eq = BloodFlowEquations2D(; h = 0.1)
xyz_data = [SA[cos(0.2*si),sin(0.2*si),si] for si in range(0,40,100)]

curve = interpolate_curve(xyz_data)
L = curve.t[end-1]
println("curve length : $L")
BloodFlowTrixi.curvature(s) = norm(DataInterpolations.derivative(curve,s,2))

mesh = P4estMesh(
    (1,2),
    polydeg = 1,
    periodicity = (true,false),
    coordinates_min = (0.0,0.0),
    coordinates_max = (2*pi,L),
    initial_refinement_level = 4
)

bc = Dict(
    :y_neg => boundary_condition_pressure_in,
    :y_pos => Trixi.BoundaryConditionDoNothing()
)

solver = DGSEM(polydeg = 1,
    surface_flux = (flux_lax_friedrichs,flux_nonconservative),
    volume_integral = VolumeIntegralFluxDifferencing((flux_lax_friedrichs,flux_nonconservative))
    )

semi = SemidiscretizationHyperbolic(
    mesh,
    eq,
    initial_condition_simple,
    source_terms = source_term_simple,
    solver,
    boundary_conditions = bc
)

tspan = (0.0, 0.3)
ode = semidiscretize(semi, tspan)

dt_adapt = StepsizeCallback(;cfl=0.5)
analyse = AliveCallback(
    alive_interval = 10,
    analysis_interval = 100
)
cb = CallbackSet(
    dt_adapt,analyse
)

sol = solve(ode, SSPRK33(),dt = dt_adapt(ode),callback= cb,saveat = 0.03,save_everystep = false)

# artery center-line

res = get3DData(eq,xyz_data,semi,sol,1)