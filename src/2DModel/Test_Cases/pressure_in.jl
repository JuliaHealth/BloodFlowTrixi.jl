function initial_condition_simple(x, t, eq::BloodFlowEquations2D; R0=2.0)
    T = eltype(x)
    A0 = T(R0^2/2)
    QRθ = T(0.0)
    Qs = T(0.0)
    E = T(1e7)
    return SVector(zero(T), QRθ,Qs, E, A0)
end

function source_term_simple(u, x, t, eq::BloodFlowEquations2D)
    T = eltype(u)
    a, Q, _, A0 = u
    A = a + A0
    s1 = zero(T)
    k = friction(u, x, eq)
    R = radius(u, eq)
    s2 = T(0.0)
    s3 = T(0.0)
    return SVector(s1, s2, s3,0, 0)
end


function boundary_condition_pressure_in(u_inner, orientation_or_normal,
direction,
x, t,
surface_flux_function,
eq::BloodFlowEquations2D)
    Pin = ifelse(t < 0.125, 2e4 * sinpi(t / 0.125)^2, 0.0)
    Ain = inv_pressure(Pin, u_inner, eq)
    A0in = u_inner[5]
    ain = Ain - A0in
    u_boundary =  SVector(
        ain,
        u_inner[2],
        u_inner[3],
        u_inner[4],
        u_inner[5]
    )
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

function boundary_condition_pressure_in(u_inner, normal,
    x, t,
    surface_flux_function,
    eq::BloodFlowEquations2D)
        Pin = ifelse(t < 0.125, 2e4 * sinpi(t / 0.125)^2, 0.0)
        A0in = u_inner[5]
        Ain = inv_pressure(Pin, u_inner, eq)
        # @info abs(Ain - A0in)
        ain = Ain - A0in
        u_boundary =  SVector(
            ain,
            u_inner[2],
            u_inner[3],
            u_inner[4],
            u_inner[5]
        )
            flux = surface_flux_function(u_boundary,u_inner, normal,
            eq)
        return flux
    end
