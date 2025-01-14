struct BloodFlowEquations1DOrd2{T <:Real,E} <: Trixi.AbstractEquationsParabolic{1, 4, GradientVariablesConservative}
    nu ::T
    model1d ::E
end
Trixi.varnames(mapin,eq::BloodFlowTrixi.BloodFlowEquations1DOrd2) = Trixi.varnames(mapin,eq.model1d)

function Trixi.flux(u,gradients,orientation::Int,eq_parab ::BloodFlowEquations1DOrd2)
    dudx = gradients
    a,Q,_,A0 = u
    A = a+A0
    val = 3*eq_parab.nu * (-(dudx[1] + dudx[4])*Q/A + dudx[2])
    return SVector(0.0,val,0,0)
end

function source_term_simple_ord2(u, x, t, eq::BloodFlowEquations1D)
    res = source_term_simple(u, x, t, eq)
    k = friction(u,x,eq)
    R = radius(u,eq)
    return SVector(res[1],res[2]/(1-R*k/(4*eq.nu)),res[3],res[4])
end

@inline function boundary_condition_pressure_in(flux_inner, u_inner,
    orientation_or_normal,direction,
    x, t,
    operator_type::Trixi.Gradient,
    equations_parabolic::BloodFlowEquations1DOrd2)
return boundary_condition_pressure_in(u_inner,orientation_or_normal,direction,x,t,flux_lax_friedrichs,equations_parabolic.model1d)
end
@inline function boundary_condition_pressure_in(flux_inner, u_inner,
    orientation_or_normal,direction,
    x, t,
    operator_type::Trixi.Divergence,
    equations_parabolic::BloodFlowEquations1DOrd2)
return flux_inner
end
# Dirichlet and Neumann boundary conditions for use with parabolic solvers in weak form.
# Note that these are general, so they apply to LaplaceDiffusion in any spatial dimension.
@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner, u_inner,
    normal::AbstractVector,
    x, t,
    operator_type::Trixi.Gradient,
    equations_parabolic::BloodFlowEquations1DOrd2)
return boundary_condition.boundary_value_function(x, t, equations_parabolic)
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner, u_inner,
    normal::AbstractVector,
    x, t,
    operator_type::Trixi.Divergence,
    equations_parabolic::BloodFlowEquations1DOrd2)
return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner, u_inner,
  normal::AbstractVector,
  x, t,
  operator_type::Trixi.Divergence,
  equations_parabolic::BloodFlowEquations1DOrd2)
return boundary_condition.boundary_normal_flux_function(x, t, equations_parabolic)
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner, u_inner,
  normal::AbstractVector,
  x, t,
  operator_type::Trixi.Gradient,
  equations_parabolic::BloodFlowEquations1DOrd2)
return flux_inner
end