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
Trixi.varnames(::typeof(cons2cons),::BloodFlowEquations1D) = ("a","Q","E","A0")

Trixi.varnames(::typeof(cons2prim),::BloodFlowEquations1D) = ("A","w","P","A0","P")

@doc raw"""
    friction(u, x, eq::BloodFlowEquations1D)

Calculates the friction term for the blood flow equations, which represents viscous resistance to flow along the artery wall.

### Parameters
- `u`: State vector containing cross-sectional area and flow rate.
- `x`: Position along the artery.
- `eq::BloodFlowEquations1D`: Instance of the blood flow model.

### Returns
Friction coefficient as a scalar.
"""
function friction(u,x,eq::BloodFlowEquations1D)
    R=radius(u,eq) 
    return eltype(u)(-11*eq.nu/R)
end

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
    pressure(u, eq::BloodFlowEquations1D)

Computes the pressure given the state vector based on the compliance of the artery.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Pressure as a scalar.
"""
function pressure(u,eq::BloodFlowEquations1D)
    T = eltype(u)
    A = u[1]+u[4]
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    b = E*h*sqrt(pi)/(1-xi^2)
    return T(b*(sqrt(A)-sqrt(A0))/A0)
end

@doc raw"""
    radius(u, eq::BloodFlowEquations1D)

Computes the radius of the artery based on the cross-sectional area.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Radius as a scalar.
"""
function radius(u,eq::BloodFlowEquations1D)
    return sqrt((u[1]+u[4])/pi)
end

@doc raw"""
    inv_pressure(p, u, eq::BloodFlowEquations1D)

Computes the inverse relation of pressure to cross-sectional area.

### Parameters
- `p`: Pressure.
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Cross-sectional area corresponding to the given pressure.
"""
function inv_pressure(p,u,eq::BloodFlowEquations1D)
    T = eltype(u)
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    # A0 p/b
    b = E*h*sqrt(pi)/(1-xi^2)
    return T((A0*p/b+sqrt(A0))^2)
end

@doc raw"""
    pressure_der(u, eq::BloodFlowEquations1D)

Computes the derivative of pressure with respect to cross-sectional area.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Derivative of pressure.
"""
function pressure_der(u,eq::BloodFlowEquations1D)
    T = eltype(u)
    A = u[1]+u[4]
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    return T(E*h*sqrt(pi)/(1-xi^2)*0.5/(sqrt(A)*A0))
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

@doc raw"""
    Trixi.cons2prim(u, eq::BloodFlowEquations1D)

Converts the conserved variables to primitive variables.

### Parameters
- `u`: State vector.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Primitive variable vector.
"""
function Trixi.cons2prim(u,eq::BloodFlowEquations1D)
    a,Q,E,A0 = u
    P = pressure(u,eq)
    A = a+A0
    w = Q/A
    return SVector(A,w,P,A0)
end

@doc raw"""
    Trixi.prim2cons(u, eq::BloodFlowEquations1D)

Converts the primitive variables to conserved variables.

### Parameters
- `u`: Primitive variable vector.
- `eq`: Instance of `BloodFlowEquations1D`.

### Returns
Conserved variable vector.
"""
function Trixi.prim2cons(u,eq::BloodFlowEquations1D)
    A,w,P,A0 = u
    a = A-A0
    Q = w*A
    E = P/sqrt(pi)*A0/(sqrt(A)-sqrt(A0))*(1-eq.xi^2)/eq.h
    return SVector(a,Q,E,A0)
end
