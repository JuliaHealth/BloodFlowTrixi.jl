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
    a_inner, Q_inner, E_inner, A0_inner = u_inner
    A_inner = a_inner+A0_inner
    PP_inner = pressure_der(u_inner, eq)
    c_inner = sqrt(A_inner*PP_inner) 
    u_equi = SVector(0.0,0.0,E_inner, A0_inner)
    PP_equi = pressure_der(u_equi, eq)
    c_equi = sqrt(A0_inner*PP_equi) 
    W2_out = Q_inner/A_inner + 4*c_inner # This should be conserve Q/A + 4 √(A P'(A))
    W1_in =  Q_inner/A_inner - 4*c_inner # This should equal the equilibrum state
    W1_out = -4*c_equi 
    A_out = inv_A_pressure_der(((W2_out - W1_out)/8)^2,u_inner,eq)
    Q_out = A_out*(W1_out+W2_out)/2
    u_boundary = SVector(A_out-A0_inner,Q_out,E_inner,A0_inner)
    # calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux1 = surface_flux_function[1](u_inner, u_boundary, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_inner, u_boundary, orientation_or_normal, eq)
    else # u_inner is "left" of boundary, u_inner is "right" of boundary
        flux1 = surface_flux_function[1](u_boundary, u_inner, orientation_or_normal, eq)
        flux2 = surface_flux_function[2](u_boundary, u_inner, orientation_or_normal, eq)
    end
    return flux1, flux2
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
