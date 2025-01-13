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

function BloodFlowEquations2D(;h,rho=1.0,xi=0.25,nu=0.04)
    return BloodFlowEquations2D(h,rho,xi,nu)
end

Trixi.have_nonconservative_terms(::BloodFlowEquations2D) = Trixi.True()

@doc raw"""
    Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations2D)
Computes the flux vector for the conservation laws of the 2D blood flow model in either the \( \theta \)-direction or the \( s \)-direction, depending on the specified orientation.

### Parameters
- `u`: State vector.
- `orientation::Integer`: Direction of the flux computation (1 for \( \theta \)-direction, otherwise \( s \)-direction).
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Flux vector as an `SVector`.
"""
function Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations2D)
    P = pressure(u, eq) # Compute pressure from state vector
    a, QRθ, Qs, E, A0 = u 
    QRθ = 0.0
    A = a + A0 # Total cross-sectional area
    if orientation == 1 # Flux in θ-direction
        f1 = QRθ / A
        f2 = QRθ^2 / (2 * A^2) + A * P
        f3 = QRθ * Qs / A^2
        return SVector(f1, f2, f3, 0, 0)
    else # Flux in s-direction
        f1 = Qs
        f2 = QRθ * Qs / A
        f3 = Qs^2 / A - QRθ^2/(2*A^2)+ A * P
        return SVector(f1, f2, f3, 0, 0)
    end
end
@doc raw"""
    Trixi.flux(u, normal, eq::BloodFlowEquations2D)
Computes the flux vector for the conservation laws of the 2D blood flow model based on a normal vector.

### Parameters
- `u`: State vector.
- `normal`: Normal vector indicating the direction of the flux.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Flux vector as an `SVector`.
"""
function Trixi.flux(u, normal, eq::BloodFlowEquations2D)
    P = pressure(u, eq) # Compute pressure from state vector
    a, QRθ, Qs, E, A0 = u 
    A = a + A0 # Total cross-sectional area
    # if normal == 1 # Flux in θ-direction
        f1 = QRθ / A
        f2 = QRθ^2 / (2 * A^2) + A * P
        f3 = QRθ * Qs / A^2
        fl1 =  SVector(f1, f2, f3, 0, 0)
    # else # Flux in s-direction
        f1 = Qs
        f2 = QRθ * Qs / A
        f3 = Qs^2 / A - QRθ^2/(2*A^2)+ A * P
        fl2 = SVector(f1, f2, f3, 0, 0)
    return fl1 .* normal[1] .+ fl2 .* normal[2]
end

@doc raw"""
    flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)
Computes the non-conservative flux for the 2D blood flow model based on the orientation.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Direction index for the flux (1 for \( \theta \)-direction, otherwise \( s \)-direction).
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

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
    flux_nonconservative(u_ll, u_rr, normal, eq::BloodFlowEquations2D)
Computes the non-conservative flux for the 2D blood flow model based on a normal vector.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `normal`: Normal vector indicating the direction of the flux.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Non-conservative flux vector.

"""
function flux_nonconservative(u_ll, u_rr, normal, eq::BloodFlowEquations2D)
    T = eltype(u_ll)
    p_ll = pressure(u_ll, eq)
    p_rr = pressure(u_rr, eq)
    pmean = (p_ll + p_rr) / 2 # Compute average pressure
    a_ll, _, _, _, A0_ll = u_ll
    a_rr, _, _, _, A0_rr = u_rr
    A_ll = a_ll + A0_ll
    A_rr = a_rr + A0_rr
    Ajump = A_rr - A_ll # Compute jump in area
    # if orientation == 1
        fn1 =  SVector(zero(T), -pmean * Ajump, 0, 0, 0)
    # else
        fn2 = SVector(zero(T), 0, -pmean * Ajump, 0, 0)
    # end
    return @. fn1*normal[1] + fn2*normal[2]
end

@doc raw"""
    Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations2D)
Computes the maximum absolute speed for wave propagation in the 2D blood flow model using a naive approach, based on the given orientation.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Direction index for the speed computation (1 for \( \theta \)-direction, otherwise \( s \)-direction).
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

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
        return max(max(abs(QRθ_ll)/A_ll^2, abs(QRθ_rr)/A_rr^2), max(sqrt(pp_ll), sqrt(pp_rr)))
    else
        ws_ll = Qs_ll / A_ll
        ws_rr = Qs_rr / A_rr
        return max(abs(ws_ll), abs(ws_rr)) + max(sqrt(A_ll * pp_ll), sqrt(A_rr * pp_rr))
    end
end

@doc raw"""
    Trixi.max_abs_speed_naive(u_ll, u_rr, normal, eq::BloodFlowEquations2D)
Computes the maximum absolute speed for wave propagation in the 2D blood flow model using a naive approach, based on a normal vector.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `normal`: Normal vector indicating the direction of wave propagation.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Maximum absolute speed.

"""
function Trixi.max_abs_speed_naive(u_ll, u_rr, normal, eq::BloodFlowEquations2D)
    a_ll, QRθ_ll, Qs_ll, _, A0_ll = u_ll
    a_rr, QRθ_rr, Qs_rr, _, A0_rr = u_rr
    A_ll = a_ll + A0_ll
    A_rr = a_rr + A0_rr
    pp_ll = pressure_der(u_ll, eq)
    pp_rr = pressure_der(u_rr, eq)
        ws_ll = Qs_ll / A_ll
        ws_rr = Qs_rr / A_rr
    return max(
        abs(ws_ll*normal[2] + sqrt(A_ll*pp_ll)*sqrt(normal[1]^2/A_ll + normal[2]^2)),
        abs(ws_rr*normal[2] + sqrt(A_rr*pp_rr)*sqrt(normal[1]^2/A_rr + normal[2]^2)),
        abs(ws_ll*normal[2] + QRθ_ll/A_ll^2*normal[1]),
        abs(ws_rr*normal[2] + QRθ_rr/A_rr^2*normal[1])
    )
end

@doc raw"""
    Trixi.max_abs_speeds(u, eq::BloodFlowEquations2D)
Computes the maximum absolute speeds for wave propagation in the 2D blood flow model in both \( \theta \)- and \( s \)-directions.

### Parameters
- `u`: State vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Tuple containing the maximum absolute speeds in the \( \theta \)- and \( s \)-directions.


"""
function Trixi.max_abs_speeds(u,eq::BloodFlowEquations2D)
    a,QRθ,Qs,E,A0 = u 
    A = a+A0
    pp= pressure_der(u,eq)
    return max(abs(QRθ/A^2),sqrt(pp)),abs(Qs/A) + sqrt(A*pp)
end

@doc raw"""
    (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations2D)
Computes the dissipation term using the Local Lax-Friedrichs method for the 2D blood flow model.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation_or_normal_direction`: Orientation index or normal vector indicating the direction of dissipation.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Dissipation vector.
"""

function (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations2D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction, eq)
    diss = -0.5 .* abs(λ) .* (u_rr .- u_ll) # Compute dissipation term
    return SVector(diss[1], diss[2], diss[3], 0, 0)
end


include("./variables.jl")
include("./bc2d.jl")
include("./Test_Cases/pressure_in.jl")
include("./Test_Cases/convergence_test.jl")