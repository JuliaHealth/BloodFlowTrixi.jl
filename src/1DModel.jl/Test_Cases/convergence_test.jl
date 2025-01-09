function initial_condition_convergence_test(x,t,eq::BloodFlowEquations1D)
    T = eltype(x)
    R0 = T(1.0)
    A0 = T(pi)*R0^2
    E = T(1e7)
    Q = T(sinpi(x[1]*t)) 
    return SVector(zero(T),Q,E,A0)
end

function source_terms_convergence_test(u,x,t,eq::BloodFlowEquations1D)
    T = eltype(u)
    A0 = u[4]
    s1 = pi*t*cospi(x[1]*t) |> T
    k = friction(u,x,eq)
    R = radius(u,eq)
    # println(2*pi*k*R*u[2]/u[1]+22*pi*eq.nu*sinpi(x[1]*t)/A0)
    s2 = pi*x[1]*cospi(x[1]*t)+ pi*t*cospi(x[1]*t)*sinpi(x[1]*t)/A0 #+2*pi*k*R*u[2]/u[1]+22*pi*eq.nu*sinpi(x[1]*t)/A0|> T
    return SVector(s1,s2,0,0)
end