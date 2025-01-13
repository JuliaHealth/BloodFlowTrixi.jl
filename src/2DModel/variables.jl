@doc raw"""
    Trixi.varnames(::typeof(cons2cons), ::BloodFlowEquations2D)

Returns the variable names in conservative form for the 2D blood flow model.

### Parameters
- `::typeof(cons2cons)`: Type representing the conservative variables.
- `::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Tuple containing the names of the conservative variables.
"""
Trixi.varnames(::typeof(cons2cons), ::BloodFlowEquations2D) = ("a", "QRθ", "Qs", "E", "A0")


@doc raw"""
    Trixi.varnames(::typeof(cons2prim), ::BloodFlowEquations2D)

Returns the variable names in primitive form for the 2D blood flow model.

### Parameters
- `::typeof(cons2prim)`: Type representing the primitive variables.
- `::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Tuple containing the names of the primitive variables.
"""
Trixi.varnames(::typeof(cons2prim), ::BloodFlowEquations2D) = ("A", "wθ", "ws", "P", "A0")


@doc raw"""
    Trixi.varnames(::typeof(cons2entropy), ::BloodFlowEquations2D)

Returns the variable names in entropy form for the 2D blood flow model.

### Parameters
- `::typeof(cons2entropy)`: Type representing the entropy variables.
- `::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Tuple containing the names of the entropy variables.
"""
Trixi.varnames(::typeof(cons2entropy), ::BloodFlowEquations2D) = ("A", "wθ", "ws", "En", "A0")


@doc raw"""
    Trixi.prim2cons(u, eq::BloodFlowEquations2D)

Converts the primitive variables to conservative variables for the 2D blood flow model.

### Parameters
- `u`: State vector in primitive form.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
State vector in conservative form as an `SVector`.
"""
function Trixi.prim2cons(u, eq::BloodFlowEquations2D)
    A, wθ, ws, P, A0 = u
    a = A - A0
    QRθ = wθ * A * sqrt(2 * A) * 3 / 4
    Qs = ws * A
    E = P * sqrt(2) * A0 / (sqrt(A) - sqrt(A0)) * (1 - eq.xi^2) / eq.h
    return SVector(a, QRθ, Qs, E, A0)
end


@doc raw"""
    Trixi.cons2prim(u, eq::BloodFlowEquations2D)

Converts the conservative variables to primitive variables for the 2D blood flow model.

### Parameters
- `u`: State vector in conservative form.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
State vector in primitive form as an `SVector`.
"""
function Trixi.cons2prim(u, eq::BloodFlowEquations2D)
    a, QRθ, Qs, E, A0 = u
    A = a + A0
    R = radius(u, eq)
    ws = Qs / A
    wθ = (4 / 3 * QRθ / R) / A
    R0 = sqrt(2 * A0)
    η = R - R0
    P = pressure(u, eq)
    return SVector(A, wθ, ws, P, A0)
end


@doc raw"""
    Trixi.cons2entropy(u, eq::BloodFlowEquations2D)

Converts the conservative variables to entropy variables for the 2D blood flow model.

### Parameters
- `u`: State vector in conservative form.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
State vector in entropy form as an `SVector`.
"""
function Trixi.cons2entropy(u, eq::BloodFlowEquations2D)
    a, QRθ, Qs, E, A0 = u
    A = a + A0
    R = radius(u, eq)
    ws = Qs / A
    wθ = (4 / 3 * QRθ / R) / A
    R0 = sqrt(2 * A0)
    η = R - R0
    P = pressure(u, eq)
    En = entropy(u, eq)
    return SVector(A, wθ, ws, En, A0)
end


@doc raw"""
    friction(u, x, eq::BloodFlowEquations2D)

Computes the friction term for the 2D blood flow model.

### Parameters
- `u`: State vector.
- `x`: Position vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Friction term as a scalar.
"""
function friction(u, x, eq::BloodFlowEquations2D)
    R = radius(u, eq) # Compute the radius based on cross-sectional area
    return eltype(u)(-11 * eq.nu / R) # Return friction term based on viscosity and radius
end


@doc raw"""
    pressure(u, eq::BloodFlowEquations2D)

Computes the pressure for the 2D blood flow model.

### Parameters
- `u`: State vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Pressure as a scalar.
"""
function pressure(u, eq::BloodFlowEquations2D)
    T = eltype(u)
    A = u[1] + u[5]
    E = u[4]
    A0 = u[5]
    R = radius(u, eq)
    R0 = sqrt(2 * A0)
    xi = eq.xi
    h = eq.h
    b = E * h / (1 - xi^2) # Precompute constant b
    return T(b * (R - R0) / R0^2)
end


@doc raw"""
    radius(u, eq::BloodFlowEquations2D)

Computes the radius based on the cross-sectional area for the 2D blood flow model.

### Parameters
- `u`: State vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Radius as a scalar.
"""
function radius(u, eq::BloodFlowEquations2D)
    return sqrt((u[1] + u[5]) * 2) # Compute radius from cross-sectional area
end

@doc raw"""
    inv_pressure(p, u, eq::BloodFlowEquations2D)

Computes the inverse of the pressure function for the 2D blood flow model.

### Parameters
- `p`: Pressure value.
- `u`: State vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Inverse pressure as a scalar.
"""
function inv_pressure(p, u, eq::BloodFlowEquations2D)
    T = eltype(u)
    E = u[4]
    A0 = u[5]
    R0 = sqrt(2 * A0)
    xi = eq.xi
    h = eq.h
    b = E * h / (1 - xi^2) # Precompute constant b
    return T((R0^2 * p / b + R0)^2 / 2)
end


@doc raw"""
    pressure_der(u, eq::BloodFlowEquations2D)

Computes the derivative of the pressure with respect to the cross-sectional area for the 2D blood flow model.

### Parameters
- `u`: State vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Derivative of pressure as a scalar.
"""
function pressure_der(u, eq::BloodFlowEquations2D)
    T = eltype(u)
    A = u[1] + u[5]
    E = u[4]
    A0 = u[5]
    xi = eq.xi
    h = eq.h
    b = E * h / (1 - xi^2)
    return T((b / sqrt(2)) * 0.5 / (sqrt(A) * A0))
end


@doc raw"""
    Trixi.entropy(u, eq::BloodFlowEquations2D)

Computes the entropy for the 2D blood flow model.

### Parameters
- `u`: State vector.
- `eq::BloodFlowEquations2D`: Instance of `BloodFlowEquations2D`.

### Returns
Entropy as a scalar.
"""
function Trixi.entropy(u, eq::BloodFlowEquations2D)
    up = cons2prim(u, eq)
    _, _, _, E, _ = u
    A, wt, ws, P, A0 = up
    psi = (ws^2 + wt^2 * 9 / 8) / 2 + P
    b = E * eq.h / (1 - eq.xi^2)
    pt = b / sqrt(2) / (3 * A0) * (A^(3 / 2) - A0^(3 / 2))
    return A * psi - pt
end
