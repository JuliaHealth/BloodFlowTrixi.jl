var documenterSearchIndex = {"docs":
[{"location":"tuto/#Tutorial-for-the-1D-model","page":"-","title":"Tutorial for the 1D model","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"In this section, we describe how to use BloodFlowTrixi.jl with Trixi.jl. This tutorial will guide you through setting up and running a 1D blood flow simulation, including mesh creation, boundary conditions, numerical fluxes, and visualization of results.","category":"page"},{"location":"tuto/#Packages","page":"-","title":"Packages","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"Before starting, ensure that the required packages are loaded:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"using Trixi\nusing BloodFlowTrixi\nusing OrdinaryDiffEq\nusing Plots","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"First, we need to choose the equation that describes the blood flow dynamics:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"eq = BloodFlowEquations1D(; h=0.1)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Here, h represents a parameter related to the initial condition or model scaling.","category":"page"},{"location":"tuto/#Mesh-and-boundary-conditions","page":"-","title":"Mesh and boundary conditions","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"We begin by defining a one-dimensional Tree mesh, which discretizes the spatial domain:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"mesh = TreeMesh(0.0, 40.0, initial_refinement_level=6, n_cells_max=10^4, periodicity=false)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This generates a non-periodic mesh for the interval 0 40, with 2^initial-refinement-level+1-1 cells. The parameter initial_refinement_level controls the initial number of cells, while n_cells_max specifies the maximum number of cells allowed during mesh refinement.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"In Trixi.jl, the Tree mesh has two labeled boundaries: ***xneg*** (left boundary) and ***xpos*** (right boundary). These labels are used to apply boundary conditions:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"bc = (\n    x_neg = boundary_condition_pressure_in,\n    x_pos = Trixi.BoundaryConditionDoNothing()\n)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"boundary_condition_pressure_in applies a pressure inflow condition at the left boundary.\nTrixi.BoundaryConditionDoNothing() specifies a \"do nothing\" boundary condition at the right boundary, meaning no flux is imposed.","category":"page"},{"location":"tuto/#Boundary-condition-implementation","page":"-","title":"Boundary condition implementation","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The inflow boundary condition is defined as:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This function applies a time-dependent pressure inflow condition.","category":"page"},{"location":"tuto/#Parameters","page":"-","title":"Parameters","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"u_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the boundary direction.\nx: Position vector.\nt: Time scalar.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.","category":"page"},{"location":"tuto/#Returns","page":"-","title":"Returns","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The boundary flux is computed based on the inflow pressure: $ P{\\text{in}} = \\begin{cases} 2 \\times 10^4 \\sin^2\\left(\\frac{\\pi t}{0.125}\\right) & \\text{if } t < 0.125 \\\n0 & \\text{otherwise} \\end{cases} $ This time-dependent inflow pressure mimics a pulsatile flow, typical in arterial blood flow. The inflow area A{\\text{in}}$ is determined using the inverse pressure relation, ensuring consistency with the physical model.","category":"page"},{"location":"tuto/#Numerical-flux","page":"-","title":"Numerical flux","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"To compute fluxes at cell interfaces, we use a combination of conservative and non-conservative fluxes:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"volume_flux = (flux_lax_friedrichs, flux_nonconservative)\nsurface_flux = (flux_lax_friedrichs, flux_nonconservative)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"flux_lax_friedrichs is a standard numerical flux for hyperbolic conservation laws.\nflux_nonconservative handles the non-conservative terms in the model, particularly those related to pressure discontinuities.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"The non-conservative flux function is defined as:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)","category":"page"},{"location":"tuto/#Parameters-2","page":"-","title":"Parameters","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"u_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations1D.","category":"page"},{"location":"tuto/#Returns-2","page":"-","title":"Returns","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The function returns the non-conservative flux vector, which is essential for capturing sharp pressure changes in the simulation.","category":"page"},{"location":"tuto/#Basis-functions-and-Shock-Capturing-DG-scheme","page":"-","title":"Basis functions and Shock Capturing DG scheme","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"To approximate the solution, we use polynomial basis functions:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"basis = LobattoLegendreBasis(2)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This defines a Lobatto-Legendre basis of polynomial degree 2, which is commonly used in high-order methods like Discontinuous Galerkin (DG) schemes.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"We then define an indicator for shock capturing, focusing on the first variable (area perturbation $a$):","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"id = IndicatorHennemannGassner(eq, basis; variable=first)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This indicator helps detect shocks or discontinuities in the solution and applies appropriate stabilization.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"The solver is defined as:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"vol = VolumeIntegralShockCapturingHG(id, volume_flux_dg=volume_flux, volume_flux_fv=surface_flux)\nsolver = DGSEM(basis, surface_flux, vol)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Here, DGSEM represents the Discontinuous Galerkin Spectral Element Method, a high-order accurate scheme suitable for hyperbolic problems.","category":"page"},{"location":"tuto/#Semi-discretization","page":"-","title":"Semi-discretization","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"We are now ready to semi-discretize the problem:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"semi = SemidiscretizationHyperbolic(\n    mesh,\n    eq,\n    initial_condition_simple,\n    source_terms = source_term_simple,\n    solver,\n    boundary_conditions=bc\n)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This step sets up the semi-discretized form of the PDE, which will be advanced in time using an ODE solver.","category":"page"},{"location":"tuto/#Source-term","page":"-","title":"Source term","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The source term accounts for additional forces acting on the blood flow, such as friction:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"source_term_simple(u, x, t, eq::BloodFlowEquations1D)","category":"page"},{"location":"tuto/#Parameters-3","page":"-","title":"Parameters","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"u: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.","category":"page"},{"location":"tuto/#Returns-3","page":"-","title":"Returns","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The source term vector is given by:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"s_1 = 0\n(no source for area perturbation).\ns_2 = frac2 pi k QR A\n, representing frictional effects.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"The friction coefficient k is computed using a model-specific friction function, and the radius R is obtained from the state vector using the radius function.","category":"page"},{"location":"tuto/#Initial-condition","page":"-","title":"Initial condition","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The initial condition specifies the starting state of the simulation:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"initial_condition_simple(x, t, eq::BloodFlowEquations1D; R0=2.0)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This function generates a simple initial condition with a uniform radius R0.","category":"page"},{"location":"tuto/#Parameters-4","page":"-","title":"Parameters","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"x: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\nR0: Initial radius (default: 2.0).","category":"page"},{"location":"tuto/#Returns-4","page":"-","title":"Returns","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The function returns a state vector with:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Zero initial area perturbation.\nZero initial flow rate.\nConstant elasticity modulus.\nReference area A_0 = pi R_0^2.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This simple initial condition is suitable for testing the model without introducing complex dynamics.","category":"page"},{"location":"tuto/#Run-the-simulation","page":"-","title":"Run the simulation","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"First, we discretize the problem in time:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Trixi.default_analysis_integrals(::BloodFlowEquations1D) = ()\ntspan = (0.0, 0.5)\node = semidiscretize(semi, tspan)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Here, tspan defines the time interval for the simulation.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Next, we add some callbacks to monitor the simulation:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"summary_callback = SummaryCallback()\nanalysis_callback = AnalysisCallback(semi, interval=200)\nstepsize_callback = StepsizeCallback(; cfl=0.5)\ncallbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"SummaryCallback provides a summary of the simulation progress.\nAnalysisCallback computes analysis metrics at specified intervals.\nStepsizeCallback adjusts the time step based on the CFL condition.","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Finally, we solve the problem:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"dt = stepsize_callback(ode)\nsol = solve(ode, SSPRK33(), dt=dt, dtmax=1e-4, dtmin=1e-11,\n            save_everystep=false, saveat=0.002, callback=callbacks)","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"Here, SSPRK33() is a third-order Strong Stability Preserving Runge-Kutta method, suitable for hyperbolic PDEs.","category":"page"},{"location":"tuto/#Plot-the-results","page":"-","title":"Plot the results","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"The results can be visualized using the following code:","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"@gif for i in eachindex(sol)\n    a1 = sol[i][1:4:end]\n    Q1 = sol[i][2:4:end]\n    A01 = sol[i][4:4:end]\n    A1 = A01 .+ a1\n    plot(Q1 ./ A1, lw=4, color=:red, ylim=(-10, 50), label=\"velocity\", legend=:bottomleft)\nend","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"This code generates an animated GIF showing the evolution of the velocity profile over time. The velocity is computed as QA, where Q is the flow rate, and A is the cross-sectional area.","category":"page"},{"location":"tuto/#Plain-code","page":"-","title":"Plain code","text":"","category":"section"},{"location":"tuto/","page":"-","title":"-","text":"using Trixi\nusing BloodFlowTrixi\nusing OrdinaryDiffEq,Plots\neq = BloodFlowEquations1D(;h=0.1)\nmesh = TreeMesh(0.0,40.0,initial_refinement_level=6,n_cells_max=10^4,periodicity=false)\nbc = (\n    x_neg = boundary_condition_pressure_in,\n    x_pos = Trixi.BoundaryConditionDoNothing()\n    )\nvolume_flux = (flux_lax_friedrichs,flux_nonconservative)\nsurface_flux = (flux_lax_friedrichs,flux_nonconservative)\nbasis = LobattoLegendreBasis(2)\nid = IndicatorHennemannGassner(eq,basis;variable=first)\nvol = VolumeIntegralShockCapturingHG(id,volume_flux_dg=volume_flux,volume_flux_fv=surface_flux)\nsolver = DGSEM(basis,surface_flux,vol)\nsemi = SemidiscretizationHyperbolic(mesh,\neq,\ninitial_condition_simple,\nsource_terms = source_term_simple,\nsolver,\nboundary_conditions=bc)\nTrixi.default_analysis_integrals(::BloodFlowEquations1D) = ()\ntspan = (0.0, 0.5)\node = semidiscretize(semi, tspan)\nsummary_callback = SummaryCallback()\nanalysis_callback = AnalysisCallback(semi, interval = 200)\nstepsize_callback = StepsizeCallback(; cfl=0.5)\ncallbacks = CallbackSet(summary_callback,analysis_callback,stepsize_callback)\ndt = stepsize_callback(ode)\nsol = solve(ode, SSPRK33(), dt = dt, dtmax = 1e-4,dtmin = 1e-11,\n            save_everystep = false,saveat = 0.002, callback = callbacks)\n\n@gif for i in eachindex(sol)\n    a1 = sol[i][1:4:end]\n    Q1 = sol[i][2:4:end]\n    A01 = sol[i][4:4:end]\n    A1 = A01.+a1\n    plot(Q1./A1,lw=4,color=:red,ylim=(-10,50),label=\"velocity\",legend=:bottomleft)\nend","category":"page"},{"location":"tuto/","page":"-","title":"-","text":"(Image: Alt Text)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = BloodFlowTrixi","category":"page"},{"location":"#BloodFlowTrixi","page":"Home","title":"BloodFlowTrixi","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BloodFlowTrixi.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [BloodFlowTrixi]","category":"page"},{"location":"#BloodFlowTrixi.BloodFlowTrixi","page":"Home","title":"BloodFlowTrixi.BloodFlowTrixi","text":"Package BloodFlowTrixi v0.0.1\n\nThis package implements 1D and 2D blood flow models for arterial circulation using Trixi.jl, enabling efficient numerical simulation and analysis.\n\nDocs under https://yolhan83.github.io/BloodFlowTrixi.jl\n\n\n\n\n\n","category":"module"},{"location":"#BloodFlowTrixi.BloodFlowEquations1D","page":"Home","title":"BloodFlowTrixi.BloodFlowEquations1D","text":"BloodFlowEquations1D(;h,rho=1.0,xi=0.25,nu=0.04)\n\nBlood Flow equations in one space dimension. This model describes the dynamics of blood flow along a compliant artery using one-dimensional equations derived from the Navier-Stokes equations. The equations account for conservation of mass and momentum, incorporating the effect of arterial compliance and frictional losses.\n\nThe governing equations are given by\n\nleftbeginaligned\n  fracpartial apartial t + fracpartialpartial x(Q) = 0 \n  fracpartial Qpartial t + fracpartialpartial xleft(fracQ^2A + A P(a)right) = P(a) fracpartial Apartial x - 2 pi R k frac Q A\n  P(a) = P_ext + fracEhsqrtpi1-xi^2fracsqrtA - sqrtA_0A_0 \n  R = sqrtfracApi\nendalignedright\n\n\n\n\n\n","category":"type"},{"location":"#BloodFlowTrixi.BloodFlowEquations2D","page":"Home","title":"BloodFlowTrixi.BloodFlowEquations2D","text":"BloodFlowEquations2D(;h,rho=1.0,xi=0.25)\n\nDefines the two-dimensional blood flow equations derived from the Navier-Stokes equations in curvilinear coordinates under the thin-artery assumption. This model describes the dynamics of blood flow along a compliant artery in two spatial dimensions (s, θ).\n\nParameters\n\nh::T: Wall thickness of the artery.\nrho::T: Fluid density (default 1.0).\nxi::T: Poisson's ratio (default 0.25).\nnu::T: Viscosity coefficient.\n\nThe governing equations account for conservation of mass and momentum, incorporating the effects of arterial compliance, curvature, and frictional losses.\n\nleftbeginaligned\n    fracpartial apartial t + fracpartialpartial thetaleft( fracQ_RthetaA right) + fracpartialpartial s(Q_s) = 0 \n    fracpartial Q_Rthetapartial t + fracpartialpartial thetaleft(fracQ_Rtheta^22A^2 + A P(a)right) + fracpartialpartial sleft( fracQ_RthetaQ_sA right) = P(a) fracpartial Apartial theta - 2 R k fracQ_RthetaA + frac2R3 mathcalCsin theta fracQ_s^2A \n    fracpartial Q_spartial t + fracpartialpartial thetaleft(fracQ_Rtheta Q_sA^2 right) + fracpartialpartial sleft( fracQ_s^2A - fracQ_Rtheta^22A^2 + A P(a) right) = P(a) fracpartial Apartial s - R k fracQ_sA - frac2R3 mathcalCsin theta fracQ_s Q_RthetaA^2 \n    P(a) = P_ext + fracEhsqrt2left(1-xi^2right)fracsqrtA - sqrtA_0A_0 \n    R = sqrt2A\nendalignedright\n\n\n\n\n\n","category":"type"},{"location":"#Trixi.DissipationLocalLaxFriedrichs-Tuple{Any, Any, Any, BloodFlowEquations1D}","page":"Home","title":"Trixi.DissipationLocalLaxFriedrichs","text":"(dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations1D)\n\nCalculates the dissipation term using the Local Lax-Friedrichs method.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation_or_normal_direction: Orientation or normal direction.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nDissipation vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.DissipationLocalLaxFriedrichs-Tuple{Any, Any, Any, BloodFlowEquations2D}","page":"Home","title":"Trixi.DissipationLocalLaxFriedrichs","text":"(dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations2D)\n\nCalculates the dissipation term using the Local Lax-Friedrichs method.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation_or_normal_direction: Orientation or normal direction.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nDissipation vector.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_outflow-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_outflow","text":"boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)\n\nImplements the outflow boundary condition, assuming that there is no reflection at the boundary.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nComputed boundary flux.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_outflow-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_outflow","text":"boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)\n\nImplements the outflow boundary condition, assuming that there is no reflection at the boundary.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nComputed boundary flux.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_pressure_in-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_pressure_in","text":"boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)\n\nImplements a pressure inflow boundary condition where the inflow pressure varies with time.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the boundary direction.\nx: Position vector.\nt: Time scalar.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nComputed boundary flux with inflow pressure specified by:\n\nP_in = begincases\n2 times 10^4 sin^2(pi t  0125)  textif  t  0125 \n0  textotherwise\nendcases\n\nThe corresponding inflow area A_{in} is computed using the inverse pressure relation, and the boundary state is constructed accordingly.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_slip_wall-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)\n\nImplements a slip wall boundary condition where the normal component of velocity is reflected.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nComputed boundary flux at the slip wall.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_slip_wall-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)\n\nImplements a slip wall boundary condition where the normal component of velocity is reflected.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nComputed boundary flux at the slip wall.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.flux_nonconservative-Tuple{Any, Any, Integer, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.flux_nonconservative","text":"flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)\n\nComputes the non-conservative flux for the model, used for handling discontinuities in pressure.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nNon-conservative flux vector.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.flux_nonconservative-Tuple{Any, Any, Integer, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.flux_nonconservative","text":"flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)\n\nComputes the non-conservative flux for the model, used for handling discontinuities in pressure.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nNon-conservative flux vector.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.friction-Tuple{Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.friction","text":"friction(u, x, eq::BloodFlowEquations1D)\n\nCalculates the friction term for the blood flow equations, which represents viscous resistance to flow along the artery wall.\n\nParameters\n\nu: State vector containing cross-sectional area and flow rate.\nx: Position along the artery.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nFriction coefficient as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.friction-Tuple{Any, Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.friction","text":"friction(u, x, eq::BloodFlowEquations2D)\n\nCalculates the friction term for the blood flow equations, representing viscous resistance to flow along the artery wall.\n\nParameters\n\nu: State vector containing cross-sectional area and flow rate.\nx: Position along the artery.\neq::BloodFlowEquations2D: Instance of the blood flow model.\n\nReturns\n\nFriction coefficient as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.initial_condition_simple-Tuple{Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.initial_condition_simple","text":"initial_condition_simple(x, t, eq::BloodFlowEquations1D; R0=2.0)\n\nGenerates a simple initial condition with a specified initial radius R0.\n\nParameters\n\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\nR0: Initial radius (default: 2.0).\n\nReturns\n\nState vector with zero initial area perturbation, zero flow rate, constant elasticity modulus, and reference area computed as A_0 = \\pi R_0^2.\n\nThis initial condition is suitable for basic tests without complex dynamics.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.inv_pressure-Tuple{Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.inv_pressure","text":"inv_pressure(p, u, eq::BloodFlowEquations1D)\n\nComputes the inverse relation of pressure to cross-sectional area.\n\nParameters\n\np: Pressure.\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nCross-sectional area corresponding to the given pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.inv_pressure-Tuple{Any, Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.inv_pressure","text":"inv_pressure(p, u, eq::BloodFlowEquations2D)\n\nComputes the inverse relation of pressure to cross-sectional area.\n\nParameters\n\np: Pressure.\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nCross-sectional area corresponding to the given pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure-Tuple{Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.pressure","text":"pressure(u, eq::BloodFlowEquations1D)\n\nComputes the pressure given the state vector based on the compliance of the artery.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nPressure as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure-Tuple{Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.pressure","text":"pressure(u, eq::BloodFlowEquations2D)\n\nComputes the pressure given the state vector based on the compliance of the artery.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nPressure as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure_der-Tuple{Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.pressure_der","text":"pressure_der(u, eq::BloodFlowEquations1D)\n\nComputes the derivative of pressure with respect to cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nDerivative of pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure_der-Tuple{Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.pressure_der","text":"pressure_der(u, eq::BloodFlowEquations2D)\n\nComputes the derivative of pressure with respect to cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nDerivative of pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.radius-Tuple{Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.radius","text":"radius(u, eq::BloodFlowEquations1D)\n\nComputes the radius of the artery based on the cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nRadius as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.radius-Tuple{Any, BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.radius","text":"radius(u, eq::BloodFlowEquations2D)\n\nComputes the radius of the artery based on the cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nRadius as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.source_term_simple-Tuple{Any, Any, Any, BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.source_term_simple","text":"source_term_simple(u, x, t, eq::BloodFlowEquations1D)\n\nComputes a simple source term for the blood flow model, focusing on frictional effects.\n\nParameters\n\nu: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nSource terms vector where:\n\ns_1 = 0 (no source for area perturbation).\ns_2 represents the friction term given by s_2 = \\frac{2 \\pi k Q}{R A}.\n\nFriction coefficient k is computed using the friction function, and the radius R is obtained using the radius function.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.cons2prim-Tuple{Any, BloodFlowEquations1D}","page":"Home","title":"Trixi.cons2prim","text":"Trixi.cons2prim(u, eq::BloodFlowEquations1D)\n\nConverts the conserved variables to primitive variables.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nPrimitive variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.cons2prim-Tuple{Any, BloodFlowEquations2D}","page":"Home","title":"Trixi.cons2prim","text":"Trixi.cons2prim(u, eq::BloodFlowEquations2D)\n\nConverts the conserved variables to primitive variables.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nPrimitive variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux-Tuple{Any, Integer, BloodFlowEquations1D}","page":"Home","title":"Trixi.flux","text":"Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations1D)\n\nComputes the flux vector for the conservation laws of the blood flow model.\n\nParameters\n\nu: State vector.\norientation::Integer: Orientation index for flux computation.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nFlux vector as an SVector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux-Tuple{Any, Integer, BloodFlowEquations2D}","page":"Home","title":"Trixi.flux","text":"Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations2D)\n\nComputes the flux vector for the conservation laws of the two-dimensional blood flow model.\n\nParameters\n\nu: State vector.\norientation::Integer: Orientation index for flux computation (1 for θ-direction, 2 for s-direction).\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nFlux vector as an SVector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.initial_condition_convergence_test-Tuple{Any, Any, BloodFlowEquations1D}","page":"Home","title":"Trixi.initial_condition_convergence_test","text":"initial_condition_convergence_test(x, t, eq::BloodFlowEquations1D)\n\nGenerates a smooth initial condition for convergence tests of the blood flow equations.\n\nParameters\n\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nInitial condition state vector with zero initial area perturbation, sinusoidal flow rate, a constant elasticity modulus, and reference area.\n\nDetails\n\nThe returned initial condition has:\n\nZero perturbation in area (a = 0).\nA sinusoidal flow rate given by Q = sin(\\pi x t).\nA constant elasticity modulus E.\nA reference cross-sectional area A_0 = \\pi R_0^2 for R_0 = 1.\n\nThis initial condition can be used to verify the accuracy and stability of numerical solvers.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.max_abs_speed_naive-Tuple{Any, Any, Integer, BloodFlowEquations1D}","page":"Home","title":"Trixi.max_abs_speed_naive","text":"Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)\n\nCalculates the maximum absolute speed for wave propagation in the blood flow model using a naive approach.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nMaximum absolute speed.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.max_abs_speed_naive-Tuple{Any, Any, Integer, BloodFlowEquations2D}","page":"Home","title":"Trixi.max_abs_speed_naive","text":"Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)\n\nCalculates the maximum absolute speed for wave propagation in the blood flow model using a naive approach.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nMaximum absolute speed.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.prim2cons-Tuple{Any, BloodFlowEquations1D}","page":"Home","title":"Trixi.prim2cons","text":"Trixi.prim2cons(u, eq::BloodFlowEquations1D)\n\nConverts the primitive variables to conserved variables.\n\nParameters\n\nu: Primitive variable vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nConserved variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.prim2cons-Tuple{Any, BloodFlowEquations2D}","page":"Home","title":"Trixi.prim2cons","text":"Trixi.prim2cons(u, eq::BloodFlowEquations2D)\n\nConverts the primitive variables to conserved variables.\n\nParameters\n\nu: Primitive variable vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nConserved variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.source_terms_convergence_test-Tuple{Any, Any, Any, BloodFlowEquations1D}","page":"Home","title":"Trixi.source_terms_convergence_test","text":"source_terms_convergence_test(u, x, t, eq::BloodFlowEquations1D)\n\nComputes the source terms for convergence tests of the blood flow equations.\n\nParameters\n\nu: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nSource terms vector.\n\nDetails\n\nThe source terms are derived based on the smooth initial condition and friction effects:\n\ns_1 represents the source term for area perturbation and is given by s_1 = \\pi t \\cos(\\pi x t).\ns_2 represents the source term for the flow rate and includes contributions from spatial and temporal variations as well as friction effects.\n\nThe radius R is computed using the radius function, and the friction coefficient k is obtained using the friction function.\n\nThis function is useful for evaluating the correctness of source term handling in numerical solvers.\n\n\n\n\n\n","category":"method"},{"location":"#Tutorial","page":"Home","title":"Tutorial","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"tuto.md\"]","category":"page"}]
}
