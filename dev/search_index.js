var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = BloodFlowTrixi","category":"page"},{"location":"#BloodFlowTrixi","page":"Home","title":"BloodFlowTrixi","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BloodFlowTrixi.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [BloodFlowTrixi]","category":"page"},{"location":"#BloodFlowTrixi.BloodFlowTrixi","page":"Home","title":"BloodFlowTrixi.BloodFlowTrixi","text":"Package BloodFlowTrixi v0.0.1\n\nThis package implements 1D and 2D blood flow models for arterial circulation using Trixi.jl, enabling efficient numerical simulation and analysis.\n\nDocs under https://yolhan83.github.io/BloodFlowTrixi.jl\n\n\n\n\n\n","category":"module"},{"location":"#BloodFlowTrixi.BloodFlowEquations1D","page":"Home","title":"BloodFlowTrixi.BloodFlowEquations1D","text":"BloodFlowEquations1D(;h,rho=1.0,xi=0.25)\n\nBlood Flow equations in one space dimension. This model describes the dynamics of blood flow along a compliant artery using one-dimensional equations derived from the Navier-Stokes equations. The equations account for conservation of mass and momentum, incorporating the effect of arterial compliance and frictional losses.\n\nThe governing equations are given by\n\nleftbeginaligned\n  fracpartial apartial t + fracpartialpartial x(Q) = 0 \n  fracpartial Qpartial t + fracpartialpartial xleft(fracQ^2A + A P(a)right) = P(a) fracpartial Apartial x - 2 pi R k frac Q A\n  P(a) = P_ext + fracEhsqrtpi1-xi^2fracsqrtA - sqrtA_0A_0 \n  R = sqrtfracApi\nendalignedright\n\n\n\n\n\n","category":"type"},{"location":"#BloodFlowTrixi.BloodFlowEquations2D","page":"Home","title":"BloodFlowTrixi.BloodFlowEquations2D","text":"BloodFlowEquations2D(;h,rho=1.0,xi=0.25)\n\nDefines the two-dimensional blood flow equations derived from the Navier-Stokes equations in curvilinear coordinates under the thin-artery assumption. This model describes the dynamics of blood flow along a compliant artery in two spatial dimensions (s, θ).\n\nParameters\n\nh::T: Wall thickness of the artery.\nrho::T: Fluid density (default 1.0).\nxi::T: Poisson's ratio (default 0.25).\nnu::T: Viscosity coefficient.\n\nThe governing equations account for conservation of mass and momentum, incorporating the effects of arterial compliance, curvature, and frictional losses.\n\nleftbeginaligned\n    fracpartial apartial t + fracpartialpartial thetaleft( fracQ_RthetaA right) + fracpartialpartial_s(Q_s) = 0 \n    fracpartial Q_Rthetapartial t + fracpartialpartial thetaleft(fracQ_Rtheta^22A^2 + A P(a)right) + fracpartialpartial_sleft( fracQ_RthetaQ_sA right) = P(a) fracpartial Apartial theta - 2 R k fracQ_RthetaA + frac2R3 mathcalCsin theta fracQ_s^2A\n    fracpartial Q_spartial t + fracpartialpartial thetaleft(fracQ_Rtheta Q_sA^2 right) + fracpartialpartial_sleft( fracQ_s^2A - fracQ_Rtheta^22A^2 + A P(a) right) = P(a) fracpartial Apartial s - R k fracQ_sA - frac2R3 mathcalCsin theta fracQ_s Q_RthetaA^2\n    P(a) = P_ext + fracEhsqrt2left(1-xi^2right)fracsqrtA - sqrtA_0A_0 \n  R = sqrt2A\nendalignedright\n\n\n\n\n\n","category":"type"},{"location":"#Trixi.DissipationLocalLaxFriedrichs-Tuple{Any, Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"Trixi.DissipationLocalLaxFriedrichs","text":"(dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations1D)\n\nCalculates the dissipation term using the Local Lax-Friedrichs method.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation_or_normal_direction: Orientation or normal direction.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nDissipation vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.DissipationLocalLaxFriedrichs-Tuple{Any, Any, Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"Trixi.DissipationLocalLaxFriedrichs","text":"(dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations2D)\n\nCalculates the dissipation term using the Local Lax-Friedrichs method.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation_or_normal_direction: Orientation or normal direction.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nDissipation vector.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_outflow-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_outflow","text":"boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)\n\nImplements the outflow boundary condition, assuming that there is no reflection at the boundary.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nComputed boundary flux.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_outflow-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_outflow","text":"boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)\n\nImplements the outflow boundary condition, assuming that there is no reflection at the boundary.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nComputed boundary flux.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_pressure_in-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_pressure_in","text":"boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)\n\nImplements a pressure inflow boundary condition where the inflow pressure varies with time.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the boundary direction.\nx: Position vector.\nt: Time scalar.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nComputed boundary flux with inflow pressure specified by:\n\nP_in = begincases\n2 times 10^4 sin^2(pi t  0125)  textif  t  0125 \n0  textotherwise\nendcases\n\nThe corresponding inflow area A_{in} is computed using the inverse pressure relation, and the boundary state is constructed accordingly.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_slip_wall-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)\n\nImplements a slip wall boundary condition where the normal component of velocity is reflected.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nComputed boundary flux at the slip wall.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.boundary_condition_slip_wall-Tuple{Any, Any, Any, Any, Any, Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)\n\nImplements a slip wall boundary condition where the normal component of velocity is reflected.\n\nParameters\n\nu_inner: State vector inside the domain near the boundary.\norientation_or_normal: Normal orientation of the boundary.\ndirection: Integer indicating the direction of the boundary.\nx: Position vector.\nt: Time.\nsurface_flux_function: Function to compute flux at the boundary.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nComputed boundary flux at the slip wall.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.flux_nonconservative-Tuple{Any, Any, Integer, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.flux_nonconservative","text":"flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)\n\nComputes the non-conservative flux for the model, used for handling discontinuities in pressure.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nNon-conservative flux vector.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.flux_nonconservative-Tuple{Any, Any, Integer, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.flux_nonconservative","text":"flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)\n\nComputes the non-conservative flux for the model, used for handling discontinuities in pressure.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nNon-conservative flux vector.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.friction-Tuple{Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.friction","text":"friction(u, x, eq::BloodFlowEquations1D)\n\nCalculates the friction term for the blood flow equations, which represents viscous resistance to flow along the artery wall.\n\nParameters\n\nu: State vector containing cross-sectional area and flow rate.\nx: Position along the artery.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nFriction coefficient as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.friction-Tuple{Any, Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.friction","text":"friction(u, x, eq::BloodFlowEquations2D)\n\nCalculates the friction term for the blood flow equations, representing viscous resistance to flow along the artery wall.\n\nParameters\n\nu: State vector containing cross-sectional area and flow rate.\nx: Position along the artery.\neq::BloodFlowEquations2D: Instance of the blood flow model.\n\nReturns\n\nFriction coefficient as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.initial_condition_convergence_test-Tuple{Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.initial_condition_convergence_test","text":"initial_condition_convergence_test(x, t, eq::BloodFlowEquations1D)\n\nGenerates a smooth initial condition for convergence tests of the blood flow equations.\n\nParameters\n\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nInitial condition state vector with zero initial area perturbation, sinusoidal flow rate, a constant elasticity modulus, and reference area.\n\nDetails\n\nThe returned initial condition has:\n\nZero perturbation in area (a = 0).\nA sinusoidal flow rate given by Q = sin(\\pi x t).\nA constant elasticity modulus E.\nA reference cross-sectional area A_0 = \\pi R_0^2 for R_0 = 1.\n\nThis initial condition can be used to verify the accuracy and stability of numerical solvers.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.initial_condition_simple-Tuple{Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.initial_condition_simple","text":"initial_condition_simple(x, t, eq::BloodFlowEquations1D; R0=2.0)\n\nGenerates a simple initial condition with a specified initial radius R0.\n\nParameters\n\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\nR0: Initial radius (default: 2.0).\n\nReturns\n\nState vector with zero initial area perturbation, zero flow rate, constant elasticity modulus, and reference area computed as A_0 = \\pi R_0^2.\n\nThis initial condition is suitable for basic tests without complex dynamics.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.inv_pressure-Tuple{Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.inv_pressure","text":"inv_pressure(p, u, eq::BloodFlowEquations1D)\n\nComputes the inverse relation of pressure to cross-sectional area.\n\nParameters\n\np: Pressure.\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nCross-sectional area corresponding to the given pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.inv_pressure-Tuple{Any, Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.inv_pressure","text":"inv_pressure(p, u, eq::BloodFlowEquations2D)\n\nComputes the inverse relation of pressure to cross-sectional area.\n\nParameters\n\np: Pressure.\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nCross-sectional area corresponding to the given pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure-Tuple{Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.pressure","text":"pressure(u, eq::BloodFlowEquations1D)\n\nComputes the pressure given the state vector based on the compliance of the artery.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nPressure as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure-Tuple{Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.pressure","text":"pressure(u, eq::BloodFlowEquations2D)\n\nComputes the pressure given the state vector based on the compliance of the artery.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nPressure as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure_der-Tuple{Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.pressure_der","text":"pressure_der(u, eq::BloodFlowEquations1D)\n\nComputes the derivative of pressure with respect to cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nDerivative of pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.pressure_der-Tuple{Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.pressure_der","text":"pressure_der(u, eq::BloodFlowEquations2D)\n\nComputes the derivative of pressure with respect to cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nDerivative of pressure.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.radius-Tuple{Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.radius","text":"radius(u, eq::BloodFlowEquations1D)\n\nComputes the radius of the artery based on the cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nRadius as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.radius-Tuple{Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"BloodFlowTrixi.radius","text":"radius(u, eq::BloodFlowEquations2D)\n\nComputes the radius of the artery based on the cross-sectional area.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nRadius as a scalar.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.source_term_simple-Tuple{Any, Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.source_term_simple","text":"source_term_simple(u, x, t, eq::BloodFlowEquations1D)\n\nComputes a simple source term for the blood flow model, focusing on frictional effects.\n\nParameters\n\nu: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nSource terms vector where:\n\ns_1 = 0 (no source for area perturbation).\ns_2 represents the friction term given by s_2 = \\frac{2 \\pi k Q}{R A}.\n\nFriction coefficient k is computed using the friction function, and the radius R is obtained using the radius function.\n\n\n\n\n\n","category":"method"},{"location":"#BloodFlowTrixi.source_terms_convergence_test-Tuple{Any, Any, Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"BloodFlowTrixi.source_terms_convergence_test","text":"source_terms_convergence_test(u, x, t, eq::BloodFlowEquations1D)\n\nComputes the source terms for convergence tests of the blood flow equations.\n\nParameters\n\nu: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.\nx: Position vector.\nt: Time scalar.\neq::BloodFlowEquations1D: Instance of the blood flow model.\n\nReturns\n\nSource terms vector.\n\nDetails\n\nThe source terms are derived based on the smooth initial condition and friction effects:\n\ns_1 represents the source term for area perturbation and is given by s_1 = \\pi t \\cos(\\pi x t).\ns_2 represents the source term for the flow rate and includes contributions from spatial and temporal variations as well as friction effects.\n\nThe radius R is computed using the radius function, and the friction coefficient k is obtained using the friction function.\n\nThis function is useful for evaluating the correctness of source term handling in numerical solvers.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.cons2prim-Tuple{Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"Trixi.cons2prim","text":"Trixi.cons2prim(u, eq::BloodFlowEquations1D)\n\nConverts the conserved variables to primitive variables.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nPrimitive variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.cons2prim-Tuple{Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"Trixi.cons2prim","text":"Trixi.cons2prim(u, eq::BloodFlowEquations2D)\n\nConverts the conserved variables to primitive variables.\n\nParameters\n\nu: State vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nPrimitive variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux-Tuple{Any, Integer, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"Trixi.flux","text":"Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations1D)\n\nComputes the flux vector for the conservation laws of the blood flow model.\n\nParameters\n\nu: State vector.\norientation::Integer: Orientation index for flux computation.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nFlux vector as an SVector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux-Tuple{Any, Integer, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"Trixi.flux","text":"Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations2D)\n\nComputes the flux vector for the conservation laws of the two-dimensional blood flow model.\n\nParameters\n\nu: State vector.\norientation::Integer: Orientation index for flux computation (1 for θ-direction, 2 for s-direction).\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nFlux vector as an SVector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.max_abs_speed_naive-Tuple{Any, Any, Integer, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"Trixi.max_abs_speed_naive","text":"Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)\n\nCalculates the maximum absolute speed for wave propagation in the blood flow model using a naive approach.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nMaximum absolute speed.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.max_abs_speed_naive-Tuple{Any, Any, Integer, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"Trixi.max_abs_speed_naive","text":"Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)\n\nCalculates the maximum absolute speed for wave propagation in the blood flow model using a naive approach.\n\nParameters\n\nu_ll: Left state vector.\nu_rr: Right state vector.\norientation::Integer: Orientation index.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nMaximum absolute speed.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.prim2cons-Tuple{Any, BloodFlowTrixi.BloodFlowEquations1D}","page":"Home","title":"Trixi.prim2cons","text":"Trixi.prim2cons(u, eq::BloodFlowEquations1D)\n\nConverts the primitive variables to conserved variables.\n\nParameters\n\nu: Primitive variable vector.\neq: Instance of BloodFlowEquations1D.\n\nReturns\n\nConserved variable vector.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.prim2cons-Tuple{Any, BloodFlowTrixi.BloodFlowEquations2D}","page":"Home","title":"Trixi.prim2cons","text":"Trixi.prim2cons(u, eq::BloodFlowEquations2D)\n\nConverts the primitive variables to conserved variables.\n\nParameters\n\nu: Primitive variable vector.\neq: Instance of BloodFlowEquations2D.\n\nReturns\n\nConserved variable vector.\n\n\n\n\n\n","category":"method"}]
}
