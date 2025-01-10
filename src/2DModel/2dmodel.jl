@doc raw"""
    BloodFlowEquations2D(;h,rho=1.0,xi=0.25)

Defines the two-dimensional blood flow equations derived from the Navier-Stokes equations in curvilinear coordinates under the thin-artery assumption. This model describes the dynamics of blood flow along a compliant artery in two spatial dimensions (s, θ).

### Parameters
- `h::T`: Wall thickness of the artery.
- `rho::T`: Fluid density (default 1.0).
- `xi::T`: Poisson's ratio (default 0.25).
- `nu::T`: Viscosity coefficient.

The governing equations account for conservation of mass and momentum, incorporating the effects of arterial compliance, curvature, and frictional losses.

```math
\left\{\begin{aligned}
    \frac{\partial a}{\partial t} + \frac{\partial}{\partial \theta}\left( \frac{Q_{R\theta}}{A} \right) + \frac{\partial}{\partial s}(Q_s) &= 0 \\
    \frac{\partial Q_{R\theta}}{\partial t} + \frac{\partial}{\partial \theta}\left(\frac{Q_{R\theta}^2}{2A^2} + A P(a)\right) + \frac{\partial}{\partial s}\left( \frac{Q_{R\theta}Q_s}{A} \right) &= P(a) \frac{\partial A}{\partial \theta} - 2 R k \frac{Q_{R\theta}}{A} + \frac{2R}{3} \mathcal{C}\sin \theta \frac{Q_s^2}{A} \\
    \frac{\partial Q_{s}}{\partial t} + \frac{\partial}{\partial \theta}\left(\frac{Q_{R\theta} Q_s}{A^2} \right) + \frac{\partial}{\partial s}\left( \frac{Q_s^2}{A} - \frac{Q_{R\theta}^2}{2A^2} + A P(a) \right) &= P(a) \frac{\partial A}{\partial s} - R k \frac{Q_s}{A} - \frac{2R}{3} \mathcal{C}\sin \theta \frac{Q_s Q_{R\theta}}{A^2} \\
    P(a) &= P_{ext} + \frac{Eh}{\sqrt{2}\left(1-\xi^2\right)}\frac{\sqrt{A} - \sqrt{A_0}}{A_0} \\
    R &= \sqrt{2A}
\end{aligned}\right.
```

"""
struct BloodFlowEquations2D{T<:Real} <: AbstractBloodFlowEquations{2,5}
    h ::T      # Wall thickness
    rho::T     # Fluid density
    xi::T      # Poisson's ratio
    nu::T      # Viscosity coefficient
end

@doc raw"""
    friction(u, x, eq::BloodFlowEquations2D)

Calculates the friction term for the blood flow equations, representing viscous resistance to flow along the artery wall.

### Parameters
- `u`: State vector containing cross-sectional area and flow rate.
- `x`: Position along the artery.
- `eq::BloodFlowEquations2D`: Instance of the blood flow model.

### Returns
Friction coefficient as a scalar.
"""
function friction(u,x,eq::BloodFlowEquations2D)
    R = radius(u,eq) # Compute the radius based on cross-sectional area
    return eltype(u)(-11 * eq.nu / R) # Return friction term based on viscosity and radius
end

@doc raw"""
    boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Implements the outflow boundary condition, assuming that there is no reflection at the boundary.

### Parameters
- `u_inner`: State vector inside the domain near the boundary.
- `orientation_or_normal`: Normal orientation of the boundary.
- `direction`: Integer indicating the direction of the boundary.
- `x`: Position vector.
- `t`: Time.
- `surface_flux_function`: Function to compute flux at the boundary.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Computed boundary flux.
"""
function boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Calculate the boundary flux without reflection
    flux = surface_flux_function(u_inner, u_inner, orientation_or_normal, eq)
    return flux
end

@doc raw"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)

Implements a slip wall boundary condition where the normal component of velocity is reflected.

### Parameters
- `u_inner`: State vector inside the domain near the boundary.
- `orientation_or_normal`: Normal orientation of the boundary.
- `direction`: Integer indicating the direction of the boundary.
- `x`: Position vector.
- `t`: Time.
- `surface_flux_function`: Function to compute flux at the boundary.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Computed boundary flux at the slip wall.
"""
function boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Create the external boundary solution state with reflected normal velocity
    u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])

    # Calculate the boundary flux based on direction
    if iseven(direction)
        flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal, eq)
    else
        flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal, eq)
    end
    return flux
end

@doc raw"""
    Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations2D)

Computes the flux vector for the conservation laws of the two-dimensional blood flow model.

### Parameters
- `u`: State vector.
- `orientation::Integer`: Orientation index for flux computation (1 for θ-direction, 2 for s-direction).
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Flux vector as an `SVector`.
"""
function Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations2D)
    P = pressure(u, eq) # Compute pressure from state vector
    a, QRθ, Qs, E, A0 = u 
    A = a + A0 # Total cross-sectional area
    if orientation == 1 # Flux in θ-direction
        f1 = QRθ / A
        f2 = QRθ^2 / (2 * A^2) + A * P
        f3 = QRθ * Qs / A^2
        return SVector(f1, f2, f3, 0, 0)
    else # Flux in s-direction
        f1 = Qs / A
        f2 = QRθ * Qs / A^2
        f3 = Qs^2 / (2 * A^2) + A * P
        return SVector(f1, f2, f3, 0, 0)
    end
end

@doc raw"""
    flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)

Computes the non-conservative flux for the model, used for handling discontinuities in pressure.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Orientation index.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Non-conservative flux vector.
"""
function flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)
    T = eltype(u_ll)
    p_ll = pressure(u_ll, eq)
    p_rr = pressure(u_rr, eq)
    pmean = (p_ll + p_rr) / 2 # Compute average pressure
    a_ll, _, _, _, A0_ll = u_ll
    a_rr, _, _, _, A0_rr = u_rr
    A_ll = a_ll + A0_ll
    A_rr = a_rr + A0_rr
    Ajump = A_rr - A_ll # Compute jump in area
    if orientation == 1
        return SVector(zero(T), -pmean * Ajump, 0, 0, 0)
    else
        return SVector(zero(T), 0, -pmean * Ajump, 0, 0)
    end
end

@doc raw"""
    Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)

Calculates the maximum absolute speed for wave propagation in the blood flow model using a naive approach.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Orientation index.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Maximum absolute speed.
"""
function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)
    a_ll, QRθ_ll, Qs_ll, _, A0_ll = u_ll
    a_rr, QRθ_rr, Qs_rr, _, A0_rr = u_rr
    A_ll = a_ll + A0_ll
    A_rr = a_rr + A0_rr
    pp_ll = pressure_der(u_ll, eq)
    pp_rr = pressure_der(u_rr, eq)
    if orientation == 1
        return max(max(abs(QRθ_ll), abs(QRθ_rr)), max(sqrt(pp_ll), sqrt(pp_rr)))
    else
        ws_ll = Qs_ll / A_ll
        ws_rr = Qs_rr / A_rr
        return max(abs(ws_ll), abs(ws_rr)) + max(sqrt(A_ll * pp_ll), sqrt(A_rr * pp_rr))
    end
end

@doc raw"""
    pressure(u, eq::BloodFlowEquations2D)

Computes the pressure given the state vector based on the compliance of the artery.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Pressure as a scalar.
"""
function pressure(u, eq::BloodFlowEquations2D)
    T = eltype(u)
    A = u[1] + u[4]
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    b = (E * h / sqrt(2)) / (1 - xi^2) # Precompute constant b
    return T(b * (sqrt(A) - sqrt(A0)) / A0)
end

@doc raw"""
    radius(u, eq::BloodFlowEquations2D)

Computes the radius of the artery based on the cross-sectional area.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Radius as a scalar.
"""
function radius(u, eq::BloodFlowEquations2D)
    return sqrt((u[1] + u[4]) / pi) # Compute radius from cross-sectional area
end

@doc raw"""
    inv_pressure(p, u, eq::BloodFlowEquations2D)

Computes the inverse relation of pressure to cross-sectional area.

### Parameters
- `p`: Pressure.
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Cross-sectional area corresponding to the given pressure.
"""
function inv_pressure(p, u, eq::BloodFlowEquations2D)
    T = eltype(u)
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    b = (E * h / sqrt(2)) / (1 - xi^2) # Precompute constant b
    return T((A0 * p / b + sqrt(A0))^2)
end

@doc raw"""
    pressure_der(u, eq::BloodFlowEquations2D)

Computes the derivative of pressure with respect to cross-sectional area.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Derivative of pressure.
"""
function pressure_der(u, eq::BloodFlowEquations2D)
    T = eltype(u)
    A = u[1] + u[4]
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    return T((E * h / sqrt(2)) / (1 - xi^2) * 0.5 / (sqrt(A) * A0))
end

@doc raw"""
    (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations2D)

Calculates the dissipation term using the Local Lax-Friedrichs method.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation_or_normal_direction`: Orientation or normal direction.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Dissipation vector.
"""
function (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations2D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction, eq)
    diss = -0.5f0 .* λ .* (u_rr .- u_ll) # Compute dissipation term
    return SVector(diss[1], diss[2], diss[3], 0, 0)
end

@doc raw"""
    Trixi.cons2prim(u, eq::BloodFlowEquations2D)

Converts the conserved variables to primitive variables.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Primitive variable vector.
"""
function Trixi.cons2prim(u, eq::BloodFlowEquations2D)
    a, QRθ, Qs, E, A0 = u
    return SVector(a, QRθ, Qs, E, A0)
end

@doc raw"""
    Trixi.prim2cons(u, eq::BloodFlowEquations2D)

Converts the primitive variables to conserved variables.

### Parameters
- `u`: Primitive variable vector.
- `eq`: Instance of `BloodFlowEquations2D`.

### Returns
Conserved variable vector.
"""
function Trixi.prim2cons(u, eq::BloodFlowEquations2D)
    a, QRθ, Qs, E, A0 = u
    return SVector(a, QRθ, Qs, E, A0)
end
