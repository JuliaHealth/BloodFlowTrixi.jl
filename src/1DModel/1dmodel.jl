@doc raw"""
    BloodFlowEquations1D(;h,rho=1.0,xi=0.25,nu=0.04)

Blood Flow equations in one space dimension. This model describes the dynamics of blood flow along a compliant artery using one-dimensional equations derived from the Navier-Stokes equations. The equations account for conservation of mass and momentum, incorporating the effect of arterial compliance and frictional losses.

The governing equations are given by
```math
\left\{\begin{aligned}
  \frac{\partial a}{\partial t} + \frac{\partial}{\partial x}(Q) &= 0 \\
  \frac{\partial Q}{\partial t} + \frac{\partial}{\partial x}\left(\frac{Q^2}{A} + A P(a)\right) &= P(a) \frac{\partial A}{\partial x} - 2 \pi R k \frac Q {A}\\
  P(a) &= P_{ext} + \frac{Eh\sqrt{\pi}}{1-\xi^2}\frac{\sqrt{A} - \sqrt{A_0}}{A_0} \\
  R &= \sqrt{\frac{A}{\pi}}
\end{aligned}\right.
```
"""
struct BloodFlowEquations1D{T<:Real} <: AbstractBloodFlowEquations{1,4}
    # constant coefficients
    h ::T      # Wall thickness
    rho::T     # Fluid density
    xi::T      # Poisson's ratio
    nu::T      # Viscosity coefficient
end

function BloodFlowEquations1D(;h,rho=1.0,xi=0.25,nu=0.04)
    return BloodFlowEquations1D(h,rho,xi,nu)
end

Trixi.have_nonconservative_terms(::BloodFlowEquations1D) = Trixi.True()

@doc raw"""
    Trixi.flux(u, orientation::Integer, eq::BloodFlowEquations1D)

Computes the flux vector for the conservation laws of the blood flow model.

### Parameters
- `u`: State vector.
- `orientation::Integer`: Orientation index for flux computation.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Flux vector as an `SVector`.
"""
function Trixi.flux(u, orientation::Integer,eq::BloodFlowEquations1D)
    # up = cons2prim(u,eq)
    P = pressure(u,eq)
    a,Q,E,A0 = u 
    A = a+A0
    f1 = Q
    f2 = Q^2/A+A*P
    return SVector(f1,f2,0,0)
end

@doc raw"""
    flux_nonconservative(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)

Computes the non-conservative flux for the model, used for handling discontinuities in pressure.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Orientation index.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Non-conservative flux vector.
"""
function flux_nonconservative(u_ll,u_rr,orientation::Integer,eq::BloodFlowEquations1D)
    T = eltype(u_ll)
    p_ll = pressure(u_ll,eq)
    p_rr = pressure(u_rr,eq)
    pmean = (p_ll+p_rr)/2
    a_ll,_,_,A0_ll = u_ll
    a_rr,_,_,A0_rr = u_rr
    A_ll = a_ll + A0_ll
    A_rr = a_rr + A0_rr
    Ajump = A_rr - A_ll
    return SVector(zero(T),-pmean*Ajump,0,0)
end

@doc raw"""
    Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, eq::BloodFlowEquations1D)

Calculates the maximum absolute speed for wave propagation in the blood flow model using a naive approach.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation::Integer`: Orientation index.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Maximum absolute speed.
"""
function Trixi.max_abs_speed_naive(u_ll,u_rr,orientation::Integer,eq ::BloodFlowEquations1D)
    a_ll,Q_ll,E_ll,A0_ll = u_ll
    a_rr,Q_rr,E_rr,A0_rr = u_rr
    A_ll = a_ll + A0_ll
    A_rr = a_rr + A0_rr
    pp_ll = pressure_der(u_ll,eq)
    pp_rr = pressure_der(u_rr,eq)
    w_ll = Q_ll/A_ll
    w_rr = Q_rr/A_rr
    return max(abs(w_ll),abs(w_rr))+max(sqrt(A_ll*pp_ll),sqrt(A_rr*pp_rr))
end

@doc raw"""
TODO
"""

function Trixi.max_abs_speeds(u,eq::BloodFlowEquations1D)
    a,Q,E,A0 = u 
    A = a+A0
    pp= pressure_der(u,eq)
    return abs(Q/A) + sqrt(A*pp)
end


@doc raw"""
    (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr, orientation_or_normal_direction, eq::BloodFlowEquations1D)

Calculates the dissipation term using the Local Lax-Friedrichs method.

### Parameters
- `u_ll`: Left state vector.
- `u_rr`: Right state vector.
- `orientation_or_normal_direction`: Orientation or normal direction.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Dissipation vector.
"""
function (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr,
    orientation_or_normal_direction,
    eq::BloodFlowEquations1D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
    eq)
    diss = -0.5f0 .* λ .* (u_rr .- u_ll)
    return SVector(diss[1], diss[2],0,0)
end

include("./1DModel/variables.jl")
include("./1DModel/bc1d.jl")
include("./1Dmodel/Test_Cases/pressure_in.jl")
include("./1Dmodel/Test_Cases/convergence_test.jl")