@doc raw"""
    initial_condition_simple(x, t, eq::BloodFlowEquations2D; R0=2.0)

Generates a simple initial condition with a specified initial radius `R0`.

### Parameters
- `x`: Position vector.
- `t`: Time scalar.
- `eq::BloodFlowEquations2D`: Instance of the blood flow model.
- `R0`: Initial radius (default: `2.0`).

### Returns
State vector with zero initial area perturbation, zero flow rate, constant elasticity modulus, and reference area computed as `A_0 = \pi R_0^2`.

This initial condition is suitable for basic tests without complex dynamics.
"""
function initial_condition_simple(x, t, eq::BloodFlowEquations2D; R0=2.0)
    T = eltype(x)
    A0 = T(R0^2/2)
    QRθ = T(0.0)
    Qs = T(0.0)
    E = T(1e7)
    return SVector(zero(T), QRθ,Qs, E, A0)
end

@doc raw"""
    source_term_simple(u, x, t, eq::BloodFlowEquations2D)

Computes a simple source term for the blood flow model, focusing on frictional effects.

### Parameters
- `u`: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.
- `x`: Position vector.
- `t`: Time scalar.
- `eq::BloodFlowEquations2D`: Instance of the blood flow model.

### Returns
Source terms vector where:
- `s_1 = 0` (no source for area perturbation).
- `s_2` represents the friction term given by `s_2 = \frac{2 \pi k Q}{R A}`.

Friction coefficient `k` is computed using the `friction` function, and the radius `R` is obtained using the `radius` function.
"""
function source_term_simple(u, x, t, eq::BloodFlowEquations2D)
    T = eltype(u)
    a, Q, _, A0 = u
    A = a + A0
    s1 = zero(T)
    k = friction(u, x, eq)
    R = radius(u, eq)
    s2 = T(0.0)
    s3 = T(0.0)
    return SVector(s1, s2, s3,0, 0)
end

@doc raw"""
    boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Implements a pressure inflow boundary condition where the inflow pressure varies with time.

### Parameters
- `u_inner`: State vector inside the domain near the boundary.
- `orientation_or_normal`: Normal orientation of the boundary.
- `direction`: Integer indicating the boundary direction.
- `x`: Position vector.
- `t`: Time scalar.
- `surface_flux_function`: Function to compute flux at the boundary.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Computed boundary flux with inflow pressure specified by:
```math
P_{in} = \begin{cases}
2 \times 10^4 \sin^2(\pi t / 0.125) & \text{if } t < 0.125 \\
0 & \text{otherwise}
\end{cases}
```
The corresponding inflow area `A_{in}` is computed using the inverse pressure relation, and the boundary state is constructed accordingly.
"""
function boundary_condition_pressure_in(u_inner, orientation_or_normal,
direction,
x, t,
surface_flux_function,
eq::BloodFlowEquations2D)
    Pin = ifelse(t < 0.125, 2e4 * sinpi(t / 0.125)^2, 0.0)
    Ain = inv_pressure(Pin, u_inner, eq)
    A0in = u_inner[5]
    ain = Ain - A0in
    u_boundary =  SVector(
        ain,
        u_inner[2],
        u_inner[3],
        u_inner[4],
        u_inner[5]
    )
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

function boundary_condition_pressure_in(u_inner, normal,
    x, t,
    surface_flux_function,
    eq::BloodFlowEquations2D)
        Pin = ifelse(t < 0.125, 2e4 * sinpi(t / 0.125)^2, 0.0)
        Ain = inv_pressure(Pin, u_inner, eq)
        A0in = u_inner[5]
        ain = Ain - A0in
        u_boundary =  SVector(
            ain,
            u_inner[2],
            u_inner[3],
            u_inner[4],
            u_inner[5]
        )
        # calculate the boundary flux
        if normal[2]>0 # u_inner is "left" of boundary, u_boundary is "right" of boundary
            flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal,
            eq)
        else # u_boundary is "left" of boundary, u_inner is "right" of boundary
            flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal,
        eq)
        end
    
        return flux
    end
