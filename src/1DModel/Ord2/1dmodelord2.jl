struct BloodFlowEquations1DOrd2{T <:Real,E} <: Trixi.AbstractEquationsParabolic{1, 2, GradientVariablesConservative}
    nu ::T
    model1d ::E
end

function Trixi.flux(u,gradients,orientation::Int,eq_parab ::BloodFlowEquations1DOrd2)
    dudx = gradients
    a,Q,_,_,A0 = u
    A = a+A0
    val = 3*eq_parab.nu * (-(dudx[1] + dudx[5])*Q/A + dudx[2])
    return SVector(0.0,val,0,0,0)
end

function boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1DOrd2)
    flux = boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq.model1d)
    return flux
end

function initial_condition_simple(x, t, eq::BloodFlowEquations1DOrd2; R0=2.0)
    return initial_condition_simple(x,t,eq.model1d;R0=R0)
end

function source_term_simple(u, x, t, eq::BloodFlowEquations1DOrd2)
    res = source_term_simple(u, x, t, eq.model1d)
    k = friction(u,x,eq.model1d)
    R = radius(u,eq.model1d)
    return SVector(res[1],res[2]/(1-R*k/4*eq.nu),res[3],res[4],res[5])
end
# Dirichlet and Neumann boundary conditions for use with parabolic solvers in weak form.
# Note that these are general, so they apply to LaplaceDiffusion in any spatial dimension.
@inline function (boundary_condition::BoundaryConditionDirichlet)(flux_inner, u_inner,
    normal::AbstractVector,
    x, t,
    operator_type::Gradient,
    equations_parabolic::BloodFlowEquations1DOrd2)
return boundary_condition.boundary_value_function(x, t, equations_parabolic)
end

@inline function (boundary_condition::BoundaryConditionDirichlet)(flux_inner, u_inner,
    normal::AbstractVector,
    x, t,
    operator_type::Divergence,
    equations_parabolic::BloodFlowEquations1DOrd2)
return flux_inner
end

@inline function (boundary_condition::BoundaryConditionNeumann)(flux_inner, u_inner,
  normal::AbstractVector,
  x, t,
  operator_type::Divergence,
  equations_parabolic::BloodFlowEquations1DOrd2)
return boundary_condition.boundary_normal_flux_function(x, t, equations_parabolic)
end

@inline function (boundary_condition::BoundaryConditionNeumann)(flux_inner, u_inner,
  normal::AbstractVector,
  x, t,
  operator_type::Gradient,
  equations_parabolic::BloodFlowEquations1DOrd2)
return flux_inner
end