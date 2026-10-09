@doc raw"""
    boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)

Implements the outflow boundary condition, assuming that there is no reflection at the boundary.

### Parameters
- `u_inner`: State vector inside the domain near the boundary.
- `orientation_or_normal`: Normal orientation of the boundary.
- `direction`: Integer indicating the direction of the boundary.
- `x`: Position vector.
- `t`: Time.
- `surface_flux_function`: Function to compute flux at the boundary.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Computed boundary flux.
"""
function boundary_condition_outflow(
    u_inner,
    orientation_or_normal,
    direction,
    x,
    t,
    surface_flux_function,
    eq::BloodFlowEquations1D,
)
    side = iseven(direction) ? 1 : -1
    u_boundary = boundary_state_outflow(u_inner, side, eq)
    # calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux1 = surface_flux_function[1](u_inner, u_boundary, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_inner, u_boundary, orientation_or_normal, eq)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux1 = surface_flux_function[1](u_boundary, u_inner, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_boundary, u_inner, orientation_or_normal, eq)
    end
    return flux1, flux2
end

@inline function boundary_state_outflow(u_inner, side, eq::BloodFlowEquations1D)
    a_inner, Q_inner, E_inner, A0_inner = u_inner
    A_inner = a_inner + A0_inner
    c_inner = sqrt(A_inner * pressure_der(u_inner, eq) / eq.rho)
    u_equilibrium = SVector(zero(a_inner), zero(Q_inner), E_inner, A0_inner)
    c_equilibrium = sqrt(A0_inner * pressure_der(u_equilibrium, eq) / eq.rho)
    # Preserve the outgoing characteristic and impose equilibrium on the incoming one.
    W_outgoing = Q_inner / A_inner + side * 4 * c_inner
    W_incoming = -side * 4 * c_equilibrium
    A_out = inv_A_pressure_der(eq.rho * ((W_outgoing - W_incoming) / 8)^2, u_inner, eq)
    Q_out = A_out * (W_outgoing + W_incoming) / 2
    return SVector(A_out - A0_inner, Q_out, E_inner, A0_inner)
end

@doc raw"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations1D)

Implements a slip wall boundary condition where the normal component of velocity is reflected.

### Parameters
- `u_inner`: State vector inside the domain near the boundary.
- `orientation_or_normal`: Normal orientation of the boundary.
- `direction`: Integer indicating the direction of the boundary.
- `x`: Position vector.
- `t`: Time.
- `surface_flux_function`: Function to compute flux at the boundary.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Computed boundary flux at the slip wall.
"""
function boundary_condition_slip_wall(
    u_inner,
    orientation_or_normal,
    direction,
    x,
    t,
    surface_flux_function,
    eq::BloodFlowEquations1D,
)
    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])

    # calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux1 = surface_flux_function[1](u_inner, u_boundary, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_inner, u_boundary, orientation_or_normal, eq)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux1 = surface_flux_function[1](u_boundary, u_inner, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_boundary, u_inner, orientation_or_normal, eq)
    end

    return flux1, flux2
end
