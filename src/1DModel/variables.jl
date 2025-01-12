Trixi.varnames(::typeof(cons2cons),::BloodFlowEquations1D) = ("a","Q","E","A0")

Trixi.varnames(::typeof(cons2prim),::BloodFlowEquations1D) = ("A","w","P","A0","P")

Trixi.varnames(::typeof(cons2entropy),::BloodFlowEquations1D) = ("A","w","En","A0","P")

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

function Trixi.cons2entropy(u,eq::BloodFlowEquations1D)
    a,Q,E,A0 = u
    P = pressure(u,eq)
    A = a+A0
    w = Q/A
    En = entropy(u,eq)
    return SVector(A,w,En,A0)
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

function Trixi.entropy(u,eq::BloodFlowEquations1D)
    up = cons2prim(u,eq)
    _,_,E,_ = u
    A,w,P,A0 = up
    psi = w^2/2+P
    b = E*eq.h/(1-eq.xi^2)
    pt = b*sqrt(pi)/(3*A0)*(A^(3/2)-A0^(3/2))
    return A*psi - pt
end