@doc raw"""
    BloodFlowEquations1D(;h,rho=1.0,xi=0.25)

Blood Flow equations in one space dimension. The equations are given by
```math
\begin{aligned}
  \frac{\partial A}{\partial t} + \frac{\partial}{\partial x}(Q) &= 0 \\
    \frac{\partial Q}{\partial t} + \frac{\partial}{\partial x}\left(\frac{Q^2}{A} + A P(A)\right) &= P(A) \frac{\partial A}{\partial x} - 2 \pi \sqrt{\frac{A}{\pi}} k \frac Q A\\
    P(A) &= P_{ext} + \frac{Eh}{1-\xi^2}\sqrt{\pi}\frac{\sqrt A - \sqrt{A_0}}{A_0} 
\end{aligned}
```
"""
struct BloodFlowEquations1D{T<:Real} <: AbstractBloodFlowEquations{1,4}
    # constant coefs
    h ::T
    rho::T
    xi::T
    nu::T
end

function BloodFlowEquations1D(;h,rho=1.0,xi=0.25,nu = 0.04)
    BloodFlowEquations1D(h,rho,xi,nu)
end

Trixi.have_nonconservative_terms(::BloodFlowEquations1D) = Trixi.True()
Trixi.varnames(::typeof(cons2cons),::BloodFlowEquations1D) = ("a","Q","E","A0")

Trixi.varnames(::typeof(cons2prim),::BloodFlowEquations1D) = ("a","Q","E","A0")


function friction(u,x,eq::BloodFlowEquations1D)
    R=radius(u,eq) 
    return eltype(u)(-11*eq.nu/R)
end




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

function Trixi.flux(u, orientation::Integer,eq::BloodFlowEquations1D)
    # up = cons2prim(u,eq)
    P = pressure(u,eq)
    a,Q,E,A0 = u 
    A = a+A0
    f1 = Q
    f2 = Q^2/A+A*P
    return SVector(f1,f2,0,0)
end

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

function Trixi.max_abs_speeds(u,eq::BloodFlowEquations1D)
    a,Q,E,A0 = u
    A = a+A0
    pp = pressure_der(u,eq)
    w = Q/A
    return (abs(w)+sqrt(A*pp),)
end

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

function radius(u,eq::BloodFlowEquations1D)
    return sqrt((u[1]+u[4])/pi)
end

function inv_pressure(p,u,eq::BloodFlowEquations1D)
    T = eltype(u)
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    # A0 p/b
    b = E*h*sqrt(pi)/(1-xi^2)
    return T((A0*p/b+sqrt(A0))^2)
    # T(b*(sqrt(A)-sqrt(A0))/A0)
end

function pressure_der(u,eq::BloodFlowEquations1D)
    T = eltype(u)
    A = u[1]+u[4]
    E = u[3]
    A0 = u[4]
    xi = eq.xi
    h = eq.h
    return T(E*h*sqrt(pi)/(1-xi^2)*0.5/(sqrt(A)*A0))
end

function (dissipation::Trixi.DissipationLocalLaxFriedrichs)(u_ll, u_rr,
    orientation_or_normal_direction,
    eq::BloodFlowEquations1D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
    eq)
    diss = -0.5f0 .* λ .* (u_rr .- u_ll)
    return SVector(diss[1], diss[2],0,0)
end

function Trixi.cons2prim(u,eq::BloodFlowEquations1D)
    a,Q,E,A0 = u
    return SVector(a,Q,E,A0)
end

function Trixi.prim2cons(u,eq::BloodFlowEquations1D)
    a,Q,E,A0 = u
    return SVector(a,Q,E,A0)
end

