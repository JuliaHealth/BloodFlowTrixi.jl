function initial_condition_simple(x,t,eq::BloodFlowEquations1D;R0=2.0)
    T = eltype(x)
    A0 = T(pi*R0^2)
    Q = T(0.0)
    E = T(1e7)
    return SVector(zero(T),Q,E,A0)
end

function source_term_simple(u,x,t,eq::BloodFlowEquations1D)
        T = eltype(u)
        a,Q,_,A0 = u
        A = a+A0 
        s1 = zero(T)
        k = friction(u,x,eq)
        R = radius(u,eq)
        s2 = T(2*pi*k/R*Q/A)
        return SVector(s1,s2,0,0)    
end

function boundary_condition_pressure_in(u_inner,orientation_or_normal,
direction,
x,t,
surface_flux_function,
eq::BloodFlowEquations1D)
    Pin = ifelse(t<0.125,2e4*sinpi(t/0.125)^2,0.0)
    Ain = inv_pressure(Pin,u_inner,eq)
    A0in = u_inner[4]
    ain = Ain - A0in
    u_boundary =  SVector(
        ain,
        u_inner[2],
        u_inner[3],
        u_inner[4]
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