


function Trixi.initial_condition_convergence_test(x, t, eq::BloodFlowEquations2D)
    T = eltype(x)
    R0 = T(1.0)
    A0 = T(R0^2/2)
    E = T(1e7)
    QRθ = 
    Qs = T(sinpi(x[1] * t))
    return SVector(zero(T), QRθ,Qs, E, A0)
end


function Trixi.source_terms_convergence_test(u, x, t, eq::BloodFlowEquations2D)
    T = eltype(u)
    A0 = u[4]
    s1 = pi * t * cospi(x[1] * t) |> T
    # k = friction(u, x, eq)
    # R = radius(u, eq)
    s2 = pi * x[1] * cospi(x[1] * t) + pi * t * cospi(x[1] * t) * sinpi(x[1] * t) / A0
    return SVector(s1, s2, 0, 0)
end
