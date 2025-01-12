Trixi.varnames(::typeof(cons2cons),::BloodFlowEquations2D) = ("a","QRθ","Qs","E","A0")

Trixi.varnames(::typeof(cons2prim),::BloodFlowEquations2D) = ("A","wtheta","ws","P","A0")


function Trixi.prim2cons(u, eq::BloodFlowEquations2D)
    A, wθ, ws, P, A0 = u
    a = A - A0
    QRθ = wθ * A * sqrt(2*A)*3/4
    Qs = ws * A
    E = P*sqrt(2)*A0/(sqrt(A)-sqrt(A0))*(1-eq.xi^2)/eq.h
    return SVector(a, QRθ, Qs, E, A0)
end


function Trixi.cons2prim(u, eq::BloodFlowEquations2D)
    a, QRθ, Qs, E, A0 = u
    A = a + A0
    R = radius(u,eq)
    ws = Qs / A
    wθ = (4/3*QRθ/R) / A
    R0 = sqrt(2*A0)
    η = R - R0
    P = pressure(u,eq)
    return SVector(A, wθ, ws, P, A0)
end

function friction(u,x,eq::BloodFlowEquations2D)
    R = radius(u,eq) # Compute the radius based on cross-sectional area
    return eltype(u)(-11 * eq.nu / R) # Return friction term based on viscosity and radius
end

function pressure(u, eq::BloodFlowEquations2D)
    T = eltype(u)
    A = u[1] + u[5]
    E = u[4]
    A0 = u[5]
    R = radius(u,eq)
    R0 = sqrt(2*A0)
    xi = eq.xi
    h = eq.h
    b = E * h / (1 - xi^2) # Precompute constant b
    return T(b * (R - R0) / R0^2)
end


function radius(u, eq::BloodFlowEquations2D)
    return sqrt((u[1] + u[5]) * 2) # Compute radius from cross-sectional area
end


function inv_pressure(p, u, eq::BloodFlowEquations2D)
    T = eltype(u)
    E = u[4]
    A0 = u[5]
    R0 = sqrt(2*A0)
    xi = eq.xi
    h = eq.h
    b = E * h / (1 - xi^2) # Precompute constant b
    return T((R0^2 * p / b + R0)^2/2)
end

function pressure_der(u, eq::BloodFlowEquations2D)
    T = eltype(u)
    A = u[1] + u[5]
    E = u[4]
    A0 = u[5]
    xi = eq.xi
    h = eq.h
    b = E*h/(1-xi^2)
    return T( (b / sqrt(2))* 0.5 / (sqrt(A) * A0))
end


function Trixi.entropy(u,eq::BloodFlowEquations2D)
    up = cons2prim(u,eq)
    _,_,_,E,_ = u
    A,wt,ws,P,A0 = up
    psi = (ws^2+wt^2*9/8)/2 + P
    b = E*eq.h/(1-eq.xi^2)
    pt = b/sqrt(2)/(3*A0)*A^(3/2)
    return A*psi - pt
end