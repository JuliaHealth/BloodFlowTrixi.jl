using Trixi
using .BloodFlowTrixi
using OrdinaryDiffEq

eq = BloodFlowEquations2D(; h = 0.1)

mesh = P4estMesh(
    (2,4),
    polydeg = 2,
    periodicity = (true,false),
    coordinates_min = (0.0,0.0),
    coordinates_max = (2*pi,40.0),
    initial_refinement_level = 4
)

bc = Dict(
    :y_neg => boundary_condition_pressure_in,
    :y_pos => Trixi.BoundaryConditionDoNothing()
)

solver = DGSEM(polydeg = 2,
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

sol = solve(ode, SSPRK33(),dt = dt_adapt(ode),callback= cb)