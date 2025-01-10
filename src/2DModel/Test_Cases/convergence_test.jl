@doc raw"""
    initial_condition_convergence_test(x, t, eq::BloodFlowEquations1D)

Generates a smooth initial condition for convergence tests of the blood flow equations.

### Parameters
- `x`: Position vector.
- `t`: Time scalar.
- `eq::BloodFlowEquations1D`: Instance of the blood flow model.

### Returns
Initial condition state vector with zero initial area perturbation, sinusoidal flow rate, a constant elasticity modulus, and reference area.

### Details
The returned initial condition has:
- Zero perturbation in area (`a = 0`).
- A sinusoidal flow rate given by `Q = sin(\pi x t)`.
- A constant elasticity modulus `E`.
- A reference cross-sectional area `A_0 = \pi R_0^2` for `R_0 = 1`.

This initial condition can be used to verify the accuracy and stability of numerical solvers.
"""
function Trixi.initial_condition_convergence_test(x, t, eq::BloodFlowEquations1D)
    T = eltype(x)
    R0 = T(1.0)
    A0 = T(pi) * R0^2
    E = T(1e7)
    Q = T(sinpi(x[1] * t))
    return SVector(zero(T), Q, E, A0)
end

@doc raw"""
    source_terms_convergence_test(u, x, t, eq::BloodFlowEquations1D)

Computes the source terms for convergence tests of the blood flow equations.

### Parameters
- `u`: State vector containing area perturbation, flow rate, elasticity modulus, and reference area.
- `x`: Position vector.
- `t`: Time scalar.
- `eq::BloodFlowEquations1D`: Instance of the blood flow model.

### Returns
Source terms vector.

### Details
The source terms are derived based on the smooth initial condition and friction effects:
- `s_1` represents the source term for area perturbation and is given by `s_1 = \pi t \cos(\pi x t)`.
- `s_2` represents the source term for the flow rate and includes contributions from spatial and temporal variations as well as friction effects.

The radius `R` is computed using the `radius` function, and the friction coefficient `k` is obtained using the `friction` function.

This function is useful for evaluating the correctness of source term handling in numerical solvers.
"""
function Trixi.source_terms_convergence_test(u, x, t, eq::BloodFlowEquations1D)
    T = eltype(u)
    A0 = u[4]
    s1 = pi * t * cospi(x[1] * t) |> T
    # k = friction(u, x, eq)
    # R = radius(u, eq)
    s2 = pi * x[1] * cospi(x[1] * t) + pi * t * cospi(x[1] * t) * sinpi(x[1] * t) / A0
    return SVector(s1, s2, 0, 0)
end
