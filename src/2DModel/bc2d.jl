function boundary_condition_outflow(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Calculate the boundary flux without reflection
    flux = surface_flux_function(u_inner, u_inner, orientation_or_normal, eq)
    return flux
end


function boundary_condition_outflow(u_inner, orientation_or_normal, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Calculate the boundary flux without reflection
    flux = surface_flux_function(u_inner, u_inner, orientation_or_normal, eq)
    return flux
end


function boundary_condition_slip_wall(u_inner, orientation_or_normal, direction, x, t, surface_flux_function, eq::BloodFlowEquations2D)
    # Create the external boundary solution state with reflected normal velocity
    u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])

    # Calculate the boundary flux based on direction
    if iseven(direction)
        flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal, eq)
    else
        flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal, eq)
    end
    return flux
end
