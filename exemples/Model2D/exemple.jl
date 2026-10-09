using Trixi
using BloodFlowTrixi
using OrdinaryDiffEqSSPRK

eq = BloodFlowEquations2D(; h=0.1)

# C=1 and R0=2 define a stress test of the reduced PDE; RC exceeds the
# nonsingular tubular geometry limit of the asymptotic derivation (RC < 1).

mesh = P4estMesh(
    (1, 2);
    polydeg=1,
    periodicity=(true, false),
    coordinates_min=(0.0, 0.0),
    coordinates_max=(2*pi, 40.0),
    initial_refinement_level=4,
)

bc = (; y_neg=boundary_condition_pressure_in, y_pos=boundary_condition_outflow)

volume_flux = (flux_central, flux_nonconservative)
surface_flux = (flux_lax_friedrichs, flux_nonconservative)
basis = LobattoLegendreBasis(1)
indicator = IndicatorHennemannGassner(eq, basis; variable=(u, eq) -> u[1] + u[5])
# Curvature 1 creates steep fronts: retain centered DG fluxes with local FV stabilization.
volume_integral = VolumeIntegralShockCapturingHG(
    indicator; volume_flux_dg=volume_flux, volume_flux_fv=surface_flux
)
solver = DGSEM(basis, surface_flux, volume_integral)

semi = SemidiscretizationHyperbolic(
    mesh,
    eq,
    initial_condition_simple,
    solver;
    source_terms=source_term_simple,
    boundary_conditions=bc,
)

tspan = (0.0, 0.3)
ode = semidiscretize(semi, tspan)

dt_adapt = StepsizeCallback(; cfl=0.5)
analyse = AliveCallback(; alive_interval=10, analysis_interval=100)
cb = CallbackSet(dt_adapt, analyse)

sol = solve(ode, SSPRK33(); dt=dt_adapt(ode), callback=cb)
