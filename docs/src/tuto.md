## Tutorial for the 1D model

In this section, we describe how to use **BloodFlowTrixi.jl** with **Trixi.jl**. This tutorial will guide you through setting up and running a 1D blood flow simulation, including mesh creation, boundary conditions, numerical fluxes, and visualization of results.

### Packages 
Before starting, ensure that the required packages are loaded:
```julia
using Trixi
using BloodFlowTrixi
using OrdinaryDiffEq
using Plots
```

First, we need to choose the equation that describes the blood flow dynamics:
```julia
eq = BloodFlowEquations1D(; h=0.1)
```
Here, `h` represents a parameter related to the initial condition or model scaling.

### Mesh and boundary conditions

We begin by defining a one-dimensional Tree mesh, which discretizes the spatial domain:
```julia
mesh = TreeMesh(0.0, 40.0, initial_refinement_level=6, n_cells_max=10^4, periodicity=false)
```
This generates a non-periodic mesh for the interval $[0, 40]$, with $2^{\text{initial_refinement_level}+1}-1$ cells. The parameter `initial_refinement_level` controls the initial number of cells, while `n_cells_max` specifies the maximum number of cells allowed during mesh refinement.

In **Trixi.jl**, the Tree mesh has two labeled boundaries: ***x_neg*** (left boundary) and ***x_pos*** (right boundary). These labels are used to apply boundary conditions:
```julia
bc = (
    x_neg = boundary_condition_pressure_in,
    x_pos = Trixi.BoundaryConditionDoNothing()
)
```
- `boundary_condition_pressure_in` applies a pressure inflow condition at the left boundary.
- `Trixi.BoundaryConditionDoNothing()` specifies a "do nothing" boundary condition at the right boundary, meaning no flux is imposed.

#### Boundary condition implementation
The inflow boundary condition is defined as:
```julia
boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)
```
This function applies a time-dependent pressure inflow condition.

#### Parameters
- `u_inner`: State vector inside the domain near the boundary.
- `orientation_or_normal`: Normal orientation of the boundary.
- `direction`: Integer indicating the boundary direction.
- `x`: Position vector.
- `t`: Time scalar.
- `surface_flux_function`: Function to compute flux at the boundary.
- `eq`: Instance of `BloodFlowEquations1D`.

#### Returns
The boundary flux is computed based on the inflow pressure:
$$
P_{\text{in}} = \begin{cases}
2 \times 10^4 \sin^2\left(\frac{\pi t}{0.125}\right) & \text{if } t < 0.125 \\
0 & \text{otherwise}
\end{cases}
$$
This time-dependent inflow pressure mimics a pulsatile flow, typical in arterial blood flow. The inflow area $A_{\text{in}}$ is determined using the inverse pressure relation, ensuring consistency with the physical model.

### Numerical flux

To compute fluxes at cell interfaces, we use a combination of conservative and non-conservative fluxes:
```julia
volume_flux = (flux_lax_friedrichs, flux_nonconservative)
surface_flux = (flux_lax_friedrichs, flux_nonconservative)
```
- `flux_lax_friedrichs` is a standard numerical flux for hyperbolic conservation laws.
- `flux_nonconservative` handles the non-conservative terms in the model, particularly those related to pressure discontinuities.

The non-conservative flux function is defined as:
```julia
flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)
```

#### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Orientation index.
- `eq`: Instance of `BloodFlowEquations1D`.

#### Returns
The function returns the non-conservative flux vector, which is essential for capturing sharp pressure changes in the simulation.

### Basis functions and Shock Capturing DG scheme

To approximate the solution, we use polynomial basis functions:
```julia
basis = LobattoLegendreBasis(2)
```
This defines a **Lobatto-Legendre** basis of polynomial degree $2$, which is commonly used in high-order methods like Discontinuous Galerkin (DG) schemes.

We then define an indicator for shock capturing, focusing on the first variable (area perturbation `$a$`):
```julia
id = IndicatorHennemannGassner(eq, basis; variable=first)
```
This indicator helps detect shocks or discontinuities in the solution and applies appropriate stabilization.

The solver is defined as:
```julia
vol = VolumeIntegralShockCapturingHG(id, volume_flux_dg=volume_flux, volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, vol)
```
Here, `DGSEM` represents the Discontinuous Galerkin Spectral Element Method, a high-order accurate scheme suitable for hyperbolic problems.

### Semi-discretization

We are now ready to semi-discretize the problem:
```julia
semi = SemidiscretizationHyperbolic(
    mesh,
    eq,
    initial_condition_simple,
    source_terms = source_term_simple,
    solver,
    boundary_conditions=bc
)
```
This step sets up the semi-discretized form of the PDE, which will be advanced in time using an ODE solver.

#### Source term
The source term accounts for additional forces acting on the blood flow, such as friction:
```julia
source_term_simple(u, x, t, eq::BloodFlowEquations1D)
```

#### Parameters
- `u`: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.
- `x`: Position vector.
- `t`: Time scalar.
- `eq::BloodFlowEquations1D`: Instance of the blood flow model.

#### Returns
The source term vector is given by:
- $s_1 = 0$ (no source for area perturbation).
- $s_2 = \frac{2 \pi k Q}{R A}$, representing frictional effects.

The friction coefficient $k$ is computed using a model-specific `friction` function, and the radius $R$ is obtained from the state vector using the `radius` function.

### Initial condition

The initial condition specifies the starting state of the simulation:
```julia
initial_condition_simple(x, t, eq::BloodFlowEquations1D; R0=2.0)
```
This function generates a simple initial condition with a uniform radius `R0`.

#### Parameters
- `x`: Position vector.
- `t`: Time scalar.
- `eq::BloodFlowEquations1D`: Instance of the blood flow model.
- `R0`: Initial radius (default: `2.0`).

#### Returns
The function returns a state vector with:
- Zero initial area perturbation.
- Zero initial flow rate.
- Constant elasticity modulus.
- Reference area $A_0 = \pi R_0^2$.

This simple initial condition is suitable for testing the model without introducing complex dynamics.

### Run the simulation

First, we discretize the problem in time:
```julia
Trixi.default_analysis_integrals(::BloodFlowEquations1D) = ()
tspan = (0.0, 0.5)
ode = semidiscretize(semi, tspan)
```
Here, `tspan` defines the time interval for the simulation.

Next, we add some callbacks to monitor the simulation:
```julia
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi, interval=200)
stepsize_callback = StepsizeCallback(; cfl=0.5)
callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback)
```
- `SummaryCallback` provides a summary of the simulation progress.
- `AnalysisCallback` computes analysis metrics at specified intervals.
- `StepsizeCallback` adjusts the time step based on the CFL condition.

Finally, we solve the problem:
```julia
dt = stepsize_callback(ode)
sol = solve(ode, SSPRK33(), dt=dt, dtmax=1e-4, dtmin=1e-11,
            save_everystep=false, saveat=0.002, callback=callbacks)
```
Here, `SSPRK33()` is a third-order Strong Stability Preserving Runge-Kutta method, suitable for hyperbolic PDEs.

### Plot the results

The results can be visualized using the following code:
```julia
@gif for i in eachindex(sol)
    a1 = sol[i][1:4:end]
    Q1 = sol[i][2:4:end]
    A01 = sol[i][4:4:end]
    A1 = A01 .+ a1
    plot(Q1 ./ A1, lw=4, color=:red, ylim=(-10, 50), label="velocity", legend=:bottomleft)
end
```
This code generates an animated GIF showing the evolution of the velocity profile over time. The velocity is computed as $Q/A$, where $Q$ is the flow rate, and $A$ is the cross-sectional area.

![Alt Text](./graph.gif)