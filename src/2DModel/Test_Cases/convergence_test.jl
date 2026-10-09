
function Trixi.initial_condition_convergence_test(x, t, eq::BloodFlowEquations2D)
    T = eltype(x)
    R0 = T(1.0)
    A0 = T(R0^2/2)
    E = T(1e7)
    QRθ = Qs = T(sinpi(x[1] * t))
    return SVector(zero(T), QRθ, Qs, E, A0)
end

function Trixi.source_terms_convergence_test(u, x, t, eq::BloodFlowEquations2D)
    T = eltype(u)
    A0 = u[5]
    Q = sinpi(x[1] * t)
    Qθ = pi * t * cospi(x[1] * t)
    Qt = pi * x[1] * cospi(x[1] * t)
    # Manufactured forcing replaces physical curvature and friction.
    return SVector(T(Qθ / A0), T(Qt + Q * Qθ / A0^2), T(Qt + 2 * Q * Qθ / A0^2), 0, 0)
end
