using StaticArrays
using Trixi

@testset "Characteristic outflow" begin
    # Return both states to inspect reconstruction and left/right ordering.
    states = ((u_ll, u_rr, orientation, eq) -> u_ll, (u_ll, u_rr, orientation, eq) -> u_rr)

    for T in (Float32, Float64), rho in (T(1), T(1.7))
        eq1 = BloodFlowEquations1D(T(0.1), rho, T(0.25), T(0.04))
        eq2 = BloodFlowEquations2D(T(0.1), rho, T(0.25), T(0.04))
        for (eq, A0, flow_index) in ((eq1, T(4pi), 2), (eq2, T(2), 3))
            @testset "$(typeof(eq))" begin
                for ratio in (T(0.8), T(1), T(1.2))
                    u = if eq isa BloodFlowEquations1D
                        SVector((ratio - 1) * A0, T(30), T(1e7), A0)
                    else
                        SVector((ratio - 1) * A0, T(0.3), T(30), T(1e7), A0)
                    end
                    A = u[1] + A0
                    App = A * BloodFlowTrixi.pressure_der(u, eq)
                    @test BloodFlowTrixi.inv_A_pressure_der(App, u, eq) ≈ A
                    @test BloodFlowTrixi.inv_A_pressure_der(App, u, eq) isa T

                    for side in (-1, 1)
                        direction = eq isa BloodFlowEquations1D ? 1 : 3
                        direction += side == 1
                        orientation = eq isa BloodFlowEquations1D ? 1 : 2
                        u_ll, u_rr = boundary_condition_outflow(
                            u, orientation, direction, nothing, T(0), states, eq
                        )
                        u_out = side == 1 ? u_rr : u_ll
                        @test (side == 1 ? u_ll : u_rr) == u
                        A_out = u_out[1] + A0
                        c = sqrt(App / rho)
                        c_out = sqrt(A_out * BloodFlowTrixi.pressure_der(u_out, eq) / rho)
                        u_rest = setindex(setindex(u, T(0), 1), T(0), flow_index)
                        c_rest = sqrt(A0 * BloodFlowTrixi.pressure_der(u_rest, eq) / rho)
                        @test u_out[flow_index] / A_out + side * 4 * c_out ≈
                            u[flow_index] / A + side * 4 * c
                        @test u_out[flow_index] / A_out - side * 4 * c_out ≈
                            -side * 4 * c_rest
                        @test u_out[(end - 1):end] == u[(end - 1):end]
                        if eq isa BloodFlowEquations2D
                            @test u_out[2] / A_out ≈ u[2] / A
                            normal = SVector(T(0), T(2 * side))
                            u_normal_inner, u_normal_out = boundary_condition_outflow(
                                u, normal, nothing, T(0), states, eq
                            )
                            @test u_normal_inner == u
                            @test u_normal_out ≈ u_out
                        end
                    end
                end

                u_rest = initial_condition_simple(SVector(T(0)), T(0), eq)
                orientation = eq isa BloodFlowEquations1D ? 1 : 2
                direction = eq isa BloodFlowEquations1D ? 2 : 4
                _, u_out = boundary_condition_outflow(
                    u_rest, orientation, direction, nothing, T(0), states, eq
                )
                @test u_out ≈ u_rest
            end
        end
    end

    @testset "1D order 2 boundary operators" begin
        eq = BloodFlowEquations1D(; h=0.1)
        eq_parabolic = BloodFlowEquations1DOrd2(eq)
        u = SVector(0.5, 30.0, 1e7, 4pi)
        flux_inner = SVector(0.0, 1.0, 0.0, 0.0)
        for direction in (1, 2)
            u_ll, u_rr = boundary_condition_outflow(
                u, 1, direction, nothing, 0.0, states, eq
            )
            u_out = iseven(direction) ? u_rr : u_ll
            @test boundary_condition_outflow(
                flux_inner, u, 1, direction, nothing, 0.0, Trixi.Gradient(), eq_parabolic
            ) ≈ u_out
            @test boundary_condition_outflow(
                flux_inner, u, 1, direction, nothing, 0.0, Trixi.Divergence(), eq_parabolic
            ) == flux_inner
        end
        u_in = boundary_condition_pressure_in(
            flux_inner, u, 1, 1, nothing, 0.0625, Trixi.Gradient(), eq_parabolic
        )
        @test length(u_in) == 4
        @test BloodFlowTrixi.pressure(u_in, eq) ≈ 2e4
        @test boundary_condition_pressure_in(
            flux_inner, u, 1, 1, nothing, 0.0625, Trixi.Divergence(), eq_parabolic
        ) == flux_inner
    end

    @testset "Order 2 diffusion dissipates flow energy" begin
        eq = BloodFlowEquations1D(; h=0.1)
        eq_parabolic = BloodFlowEquations1DOrd2(eq)
        mesh = TreeMesh(0.0, 2pi; initial_refinement_level=3, periodicity=true)
        solver = DGSEM(;
            polydeg=3,
            surface_flux=(flux_lax_friedrichs, flux_nonconservative),
            volume_integral=VolumeIntegralFluxDifferencing((
                flux_central, flux_nonconservative
            )),
        )
        initial_condition(x, t, eq) = SVector(0.0, sin(x[1]), 1e7, 4pi)
        semi = SemidiscretizationHyperbolicParabolic(
            mesh,
            (eq, eq_parabolic),
            initial_condition,
            solver;
            boundary_conditions=(boundary_condition_periodic, boundary_condition_periodic),
        )
        ode = semidiscretize(semi, (0.0, 0.1))
        du = similar(ode.u0)
        Trixi.rhs_parabolic!(du, ode.u0, semi, 0.0)
        u_nodes = reshape(ode.u0, 4, 4, :)
        du_nodes = reshape(du, 4, 4, :)
        energy_derivative = sum(
            solver.basis.weights[i] * u_nodes[2, i, element] * du_nodes[2, i, element] for
            element in axes(u_nodes, 3), i in axes(u_nodes, 2)
        )
        @test energy_derivative < 0
        @test all(iszero, du_nodes[[1, 3, 4], :, :])
    end

    @testset "Curvature exchanges momentum without adding kinetic energy" begin
        eq = BloodFlowEquations2D(; h=0.1, nu=0.0)
        u = SVector(0.2, 0.3, 30.0, 1e7, 2.0)
        A = u[1] + u[5]
        source = source_term_simple(u, SVector(pi / 2, 0.0), 0.0, eq)
        @test u[2] / A^2 * source[2] + u[3] / A * source[3] ≈ 0.0 atol=1e-12
    end

    @testset "2D flux orientations agree with normal fluxes" begin
        eq = BloodFlowEquations2D(; h=0.1)
        u = SVector(0.2, 0.3, 30.0, 1e7, 2.0)
        for orientation in (1, 2)
            normal = orientation == 1 ? SVector(1.0, 0.0) : SVector(0.0, 1.0)
            @test Trixi.flux(u, orientation, eq) ≈ Trixi.flux(u, normal, eq)
            speed = Trixi.max_abs_speed_naive(u, u, orientation, eq)
            @test Trixi.max_abs_speed_naive(u, u, normal, eq) ≈ speed
            @test Trixi.max_abs_speed_naive(u, u, -normal, eq) ≈ speed
            @test Trixi.max_abs_speed_naive(u, u, 2 * normal, eq) ≈ 2 * speed
        end
    end
end
