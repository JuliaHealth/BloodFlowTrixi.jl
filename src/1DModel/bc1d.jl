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
function boundary_condition_outflow(u_inner, orientation_or_normal,       
    direction,
    x, t,
    surface_flux_function,
    eq::BloodFlowEquations1D)
        # calculate the boundary flux
        if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
            flux = surface_flux_function(u_inner, u_inner, orientation_or_normal,
            eq)
            else # u_inner is "left" of boundary, u_inner is "right" of boundary
            flux = surface_flux_function(u_inner, u_inner, orientation_or_normal,
            eq)
            end
    return flux
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
function boundary_condition_slip_wall(u_inner, orientation_or_normal,       
    direction,
    x, t,
    surface_flux_function,
    eq::BloodFlowEquations1D)
    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
    -u_inner[2],
    u_inner[3],
    u_inner[4])

    # calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
    flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal,
    eq)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
    flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal,
    eq)
    end

    return flux
end