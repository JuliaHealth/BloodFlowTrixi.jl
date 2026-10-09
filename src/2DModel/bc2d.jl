@doc raw"""
    boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies a characteristic outflow condition at axial boundaries, preserving the outgoing
characteristic and prescribing the incoming characteristic from the rest state.
Circumferential boundaries use extrapolation.

### Parameters
- `u_inner`: Inner state vector at the boundary.
- `orientation_or_normal`: Orientation index or normal vector indicating the boundary direction.
- `direction`: Boundary side (1/2 for negative/positive \( \theta \), 3/4 for negative/positive \( s \)).
- `x`: Position vector at the boundary.
- `t`: Time value.
- `surface_flux_function`: Function to compute the surface flux.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Boundary flux as an `SVector`.
"""
function boundary_condition_outflow(
    u_inner,
    orientation_or_normal,
    direction,
    x,
    t,
    surface_flux_function,
    eq::BloodFlowEquations2D,
)
    side = iseven(direction) ? 1 : -1
    u_boundary =
        orientation_or_normal == 2 ? boundary_state_outflow(u_inner, side, eq) : u_inner
    if iseven(direction)
        flux1 = surface_flux_function[1](u_inner, u_boundary, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_inner, u_boundary, orientation_or_normal, eq)
    else
        flux1 = surface_flux_function[1](u_boundary, u_inner, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_boundary, u_inner, orientation_or_normal, eq)
    end
    return flux1, flux2
end

@doc raw"""
    boundary_condition_outflow(u_inner, orientation_or_normal, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies characteristic outflow at boundaries whose normal is aligned with \( s \).
Other normals use extrapolation. The normal may be scaled, as on a `P4estMesh`.

### Parameters
- `u_inner`: Inner state vector at the boundary.
- `orientation_or_normal`: Orientation index or normal vector indicating the boundary direction.
- `x`: Position vector at the boundary.
- `t`: Time value.
- `surface_flux_function`: Function to compute the surface flux.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Boundary flux as an `SVector`.
"""
function boundary_condition_outflow(
    u_inner, orientation_or_normal, x, t, surface_flux_function, eq::BloodFlowEquations2D
)
    normal = orientation_or_normal
    side = normal[2] > 0 ? 1 : -1
    u_boundary = if iszero(normal[1]) && !iszero(normal[2])
        boundary_state_outflow(u_inner, side, eq)
    else
        u_inner
    end
    flux1 = surface_flux_function[1](u_inner, u_boundary, normal, eq)
    flux2 = surface_flux_function[2](u_inner, u_boundary, normal, eq)
    return flux1, flux2
end

@inline function boundary_state_outflow(u_inner, side, eq::BloodFlowEquations2D)
    a_inner, QRθ_inner, Qs_inner, E_inner, A0_inner = u_inner
    A_inner = a_inner + A0_inner
    c_inner = sqrt(A_inner * pressure_der(u_inner, eq) / eq.rho)
    u_equilibrium = SVector(
        zero(a_inner), zero(QRθ_inner), zero(Qs_inner), E_inner, A0_inner
    )
    c_equilibrium = sqrt(A0_inner * pressure_der(u_equilibrium, eq) / eq.rho)
    W_outgoing = Qs_inner / A_inner + side * 4 * c_inner
    W_incoming = -side * 4 * c_equilibrium
    A_out = inv_A_pressure_der(eq.rho * ((W_outgoing - W_incoming) / 8)^2, u_inner, eq)
    Qs_out = A_out * (W_outgoing + W_incoming) / 2
    # QRθ / A is transported by the axial flow and remains constant across acoustic waves.
    QRθ_out = QRθ_inner * A_out / A_inner
    return SVector(A_out - A0_inner, QRθ_out, Qs_out, E_inner, A0_inner)
end

@doc raw"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies a slip-wall boundary condition for the 2D blood flow model by reflecting the normal component of the velocity at the boundary.

### Parameters
- `u_inner`: Inner state vector at the boundary.
- `orientation_or_normal`: Orientation index or normal vector indicating the boundary direction.
- `direction`: Index indicating the spatial direction (1 for \( \theta \)-direction, otherwise \( s \)-direction).
- `x`: Position vector at the boundary.
- `t`: Time value.
- `surface_flux_function`: Function to compute the surface flux.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Boundary flux as an `SVector`.
"""
function boundary_condition_slip_wall(
    u_inner,
    orientation_or_normal,
    direction,
    x,
    t,
    surface_flux_function,
    eq::BloodFlowEquations2D,
)
    # Create the external boundary solution state with reflected normal velocity
    u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])

    # Calculate the boundary flux based on direction
    if iseven(direction)
        flux1 = surface_flux_function[1](u_inner, u_boundary, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_inner, u_boundary, orientation_or_normal, eq)
    else
        flux1 = surface_flux_function[1](u_boundary, u_inner, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_boundary, u_inner, orientation_or_normal, eq)
    end
    return flux1, flux2
end
