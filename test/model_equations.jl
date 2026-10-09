using ForwardDiff
using LinearAlgebra
import BloodFlowTrixi: pressure

@testset "Published model equations" begin
    for rho in (1.0, 1.7), ratio in (0.8, 1.0, 1.2)
        eq1 = BloodFlowEquations1D(; h=0.1, rho)
        eq2 = BloodFlowEquations2D(; h=0.1, rho)
        for (eq, A0) in ((eq1, 4pi), (eq2, 2.0))
            @testset "$(typeof(eq)), rho=$rho, A/A0=$ratio" begin
                u = eq isa BloodFlowEquations1D ?
                    SVector((ratio - 1) * A0, 30.0, 1e7, A0) :
                    SVector((ratio - 1) * A0, 0.3, 30.0, 1e7, A0)
                A = u[1] + A0
                R = eq isa BloodFlowEquations1D ? sqrt(A / pi) : sqrt(2A)
                R0 = eq isa BloodFlowEquations1D ? sqrt(A0 / pi) : sqrt(2A0)
                P = u[end - 1] * eq.h / (1 - eq.xi^2) * (R - R0) / R0^2
                @test pressure(u, eq) ≈ P atol=1e-10
                @test BloodFlowTrixi.inv_pressure(P, u, eq) ≈ A
                if ratio != 1.0 # Zero pressure at rest does not determine E.
                    @test prim2cons(cons2prim(u, eq), eq) ≈ u
                end
                dp = ForwardDiff.gradient(v -> pressure(v, eq), u)[1]
                @test BloodFlowTrixi.pressure_der(u, eq) ≈ dp
                @test BloodFlowTrixi.inv_A_pressure_der(A * dp, u, eq) ≈ A
                p = P / rho
                pp = dp / rho
                orientations = eq isa BloodFlowEquations1D ? (1,) : (1, 2)
                gradient_energy = ForwardDiff.gradient(v -> entropy(v, eq), u)
                for orientation in orientations
                    # The continuous operator is dF/dU plus the infinitesimal
                    # nonconservative jump, including frozen E and A0 columns.
                    H = ForwardDiff.jacobian(v -> flux(v, orientation, eq), u) +
                        ForwardDiff.jacobian(v -> flux_nonconservative(u, v, orientation, eq), u)
                    if eq isa BloodFlowEquations1D
                        w = u[2] / A
                        expected = sort([w - sqrt(A * pp), w + sqrt(A * pp)])
                        energy_flux = v -> begin
                            Av = v[1] + v[end]
                            v[2] / Av * (v[2]^2 / (2Av) + Av * pressure(v, eq) / rho)
                        end
                    else
                        w = u[3] / A
                        expected = orientation == 1 ?
                            sort([-sqrt(pp), sqrt(pp), u[2] / A^2]) :
                            sort([w - sqrt(A * pp), w, w + sqrt(A * pp)])
                        energy_flux = v -> begin
                            Av = v[1] + v[end]
                            speed = orientation == 1 ? v[2] / Av^2 : v[3] / Av
                            speed * (v[3]^2 / (2Av) + Av * pressure(v, eq) / rho)
                        end
                    end
                    n = length(expected)
                    @test sort(eigvals(Matrix(H[1:n, 1:n]))) ≈ expected
                    @test vec(gradient_energy' * H) ≈ ForwardDiff.gradient(energy_flux, u)
                    @test max_abs_speed_naive(u, u, orientation, eq) ≈ maximum(abs, expected)
                    @test Trixi.max_abs_speeds(u, eq)[orientation] ≈ maximum(abs, expected)
                end

                source = source_term_simple(u, SVector(pi / 2, 0.0), 0.0, eq)
                k = -11 * eq.nu / R
                if eq isa BloodFlowEquations1D
                    @test dot(gradient_energy, source) ≈ 2pi * R * k * (u[2] / A)^2
                    source_ord2 = source_term_simple_ord2(u, nothing, 0.0, eq)
                    @test dot(gradient_energy, source_ord2) ≈
                        2pi * R * k / (1 - R * k / (4eq.nu)) * (u[2] / A)^2
                    gradients = (SVector(0.4, 2.0, 0.0, -0.1),)
                    diffusion = flux(u, gradients, 1, BloodFlowEquations1DOrd2(eq))
                    @test diffusion[2] ≈ 3eq.nu * (2 - (0.4 - 0.1) * u[2] / A)
                else
                    for normal in (SVector(0.7, -1.2), SVector(-2.0, 0.3))
                        Hn = ForwardDiff.jacobian(v -> flux(v, normal, eq), u) +
                            ForwardDiff.jacobian(v -> flux_nonconservative(u, v, normal, eq), u)
                        ws = u[3] / A
                        acoustic = sqrt(pp * (normal[1]^2 + A * normal[2]^2))
                        expected = sort([normal[2] * ws - acoustic,
                                         normal[2] * ws + acoustic,
                                         normal[2] * ws + normal[1] * u[2] / A^2])
                        @test sort(eigvals(Matrix(Hn[1:3, 1:3]))) ≈ expected
                        @test max_abs_speed_naive(u, u, normal, eq) >= maximum(abs, expected)
                    end
                    wt = 4u[2] / (3R * A)
                    ws = u[3] / A
                    # Curvature cancels from energy production; this coefficient
                    # detects 3Rk instead of the published 2Rk in angular friction.
                    @test dot(gradient_energy, source) ≈ 9 / 4 * R * k * wt^2 + R * k * ws^2
                    kinetic = u[2]^2 / (2A^2) + u[3]^2 / (2A)
                    beta = u[end - 1] * eq.h / (rho * (1 - eq.xi^2) * sqrt(2))
                    @test entropy(u, eq) ≈ kinetic + A * p -
                        beta / (3A0) * (A^(3 / 2) - A0^(3 / 2))
                end
            end
        end
    end

    @testset "Zero-viscosity limit" begin
        eq = BloodFlowEquations1D(; h=0.1, nu=0.0)
        u = SVector(0.2, 30.0, 1e7, 4pi)
        @test all(iszero, source_term_simple_ord2(u, nothing, 0.0, eq))
        @test all(iszero, flux(u, (ones(typeof(u)),), 1, BloodFlowEquations1DOrd2(eq)))
    end

    @testset "Manufactured source satisfies the PDE" begin
        x = SVector(0.4, 0.7)
        t = 0.3
        for eq in (BloodFlowEquations1D(; h=0.1), BloodFlowEquations2D(; h=0.1))
            exact(y, time) = initial_condition_convergence_test(y .* one(time), time, eq)
            u = exact(x, t)
            # Differentiate in time using a one-coordinate static vector.
            residual = ForwardDiff.jacobian(tau -> exact(x, tau[1]), SVector(t))[:, 1]
            for orientation in (eq isa BloodFlowEquations1D ? (1,) : (1, 2))
                jac_flux = ForwardDiff.jacobian(y -> flux(exact(y, t), orientation, eq), x)
                jac_noncons = ForwardDiff.jacobian(
                    y -> flux_nonconservative(u, exact(y, t), orientation, eq), x
                )
                residual += jac_flux[:, orientation] + jac_noncons[:, orientation]
            end
            @test source_terms_convergence_test(u, x, t, eq) ≈ residual
        end
    end

    @testset "Visualization uses physical pressure and velocity" begin
        eq1 = BloodFlowEquations1D(; h=0.1, rho=1.7)
        u1 = SVector(0.2, 30.0, 1e7, 4pi)
        semi1 = (; cache=(; elements=(; node_coordinates=[0.0])))
        data1 = get3DData(eq1, semi1, [collect(u1)]; theta_disc=3)
        @test all(P -> P ≈ pressure(u1, eq1), data1.P)
        eq2 = BloodFlowEquations2D(; h=0.1, rho=1.7)
        u2 = SVector(0.2, 0.3, 30.0, 1e7, 2.0)
        semi2 = (; cache=(; elements=(; node_coordinates=reshape([0.0, 0.0], 2, 1, 1, 1))))
        data2 = get3DData(eq2, semi2, (; u=[collect(u2)]))
        @test data2.wtheta[1] ≈ cons2prim(u2, eq2)[2]
        @test data2.ws[1] ≈ cons2prim(u2, eq2)[3]
    end
end
