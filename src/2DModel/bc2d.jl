@doc raw"""
    boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies an outflow boundary condition for the 2D blood flow model without reflecting any flux.

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
function boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Calculate the boundary flux without reflection
    flux1 = surface_flux_function[1](u_inner, u_inner, orientation_or_normal, eq)
    flux2 = surface_flux_function[2](u_inner, u_inner, orientation_or_normal, eq)
    return flux1,flux2
end


@doc raw"""
    boundary_condition_outflow(u_inner, orientation_or_normal, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies an outflow boundary condition for the 2D blood flow model without reflecting any flux. This version does not use a specific direction parameter.

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
function boundary_condition_outflow(u_inner, orientation_or_normal, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Calculate the boundary flux without reflection
    flux1 = surface_flux_function[1](u_inner, u_inner, orientation_or_normal, eq)
    flux2 = surface_flux_function[2](u_inner, u_inner, orientation_or_normal, eq)
    return flux1,flux2
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
function boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
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
    return flux1,flux2
end
