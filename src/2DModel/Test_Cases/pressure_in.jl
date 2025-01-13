@doc raw"""
    initial_condition_simple(x, t, eq::BloodFlowEquations2D; R0=2.0)

Defines a simple initial condition for the 2D blood flow model.

### Parameters
- `x`: Position vector.
- `t`: Initial time.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.
- `R0`: Initial radius (default is 2.0).

### Returns
State vector as an `SVector`.
"""
function initial_condition_simple(x, t, eq::BloodFlowEquations2D; R0=2.0)
    T = eltype(x)
    A0 = T(R0^2 / 2)
    QRθ = T(0.0)
    Qs = T(0.0)
    E = T(1e7)
    return SVector(zero(T), QRθ, Qs, E, A0)
end


@doc raw"""
    curvature(x)

Returns a constant curvature for the 2D blood flow model.

### Parameters
- `x`: Position vector.

### Returns
Curvature as a scalar.
"""
curvature(x) = typeof(x)(1.0)


@doc raw"""
    source_term_simple(u, x, t, eq::BloodFlowEquations2D)

Computes a simple source term for the 2D blood flow model, including friction and curvature effects.

### Parameters
- `u`: State vector.
- `x`: Position vector.
- `t`: Time value.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Source term as an `SVector`.
"""
function source_term_simple(u, x, t, eq::BloodFlowEquations2D)
    T = eltype(u)
    a, QRθ, Qs, _, A0 = u
    A = a + A0
    s1 = zero(T)
    k = friction(u, x, eq)
    R = radius(u, eq)
    s2 = T(
        2 * R / 3 * curvature(x[2]) * sin(x[1]) * Qs^2 / A + 3 * R * k * QRθ / A
    )
    s3 = T(
        -2 * R / 3 * curvature(x[2]) * sin(x[1]) * Qs * QRθ / A + R * k * Qs / A
    )
    return SVector(s1, s2, s3, 0, 0)
end


@doc raw"""
    boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies an inflow boundary condition with a prescribed pressure for the 2D blood flow model.

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
function boundary_condition_pressure_in(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    Pin = ifelse(t < 0.125, 2e4 * sinpi(t / 0.125)^2, 0.0)
    Ain = inv_pressure(Pin, u_inner, eq)
    A0in = u_inner[5]
    ain = Ain - A0in
    u_boundary = SVector(
        ain,
        u_inner[2],
        u_inner[3],
        u_inner[4],
        u_inner[5]
    )
    # Calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal, eq)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal, eq)
    end

    return flux
end


@doc raw"""
    boundary_condition_pressure_in(u_inner, normal, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Applies an inflow boundary condition with a prescribed pressure for the 2D blood flow model. This version does not use a specific direction parameter.

### Parameters
- `u_inner`: Inner state vector at the boundary.
- `normal`: Normal vector indicating the boundary direction.
- `x`: Position vector at the boundary.
- `t`: Time value.
- `surface_flux_function`: Function to compute the surface flux.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Boundary flux as an `SVector`.
"""
function boundary_condition_pressure_in(u_inner, normal, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    Pin = ifelse(t < 0.125, 2e4 * sinpi(t / 0.125)^2, 0.0)
    A0in = u_inner[5]
    Ain = inv_pressure(Pin, u_inner, eq)
    ain = Ain - A0in
    u_boundary = SVector(
        ain,
        u_inner[2],
        u_inner[3],
        u_inner[4],
        u_inner[5]
    )
    flux = surface_flux_function(u_inner, u_boundary, normal, eq)
    return flux
end
