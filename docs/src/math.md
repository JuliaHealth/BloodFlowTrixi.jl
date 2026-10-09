# Blood flow equations

The models follow [Mannes et al., 1D](https://doi.org/10.4236/jamp.2025.1310198)
and [Mannes et al., 2D](https://doi.org/10.4236/jamp.2025.1311220).
## Pressure and material variables

The stored state uses a = A - A0. Young's modulus E and reference area A0
are stationary fields. The function pressure returns physical transmural
pressure P, with zero external pressure as reference. Below, p = P/rho is
specific pressure, as in the derivation sections of the papers.
The function pressure_der differentiates P, and inv_A_pressure_der inverts
A * pressure_der without a density factor.

The code chooses the non-positive Navier coefficient k = -11nu/R.
This is a friction law choice; the reduced equations allow other negative k.

## One-dimensional models

Here A = pi R^2, Q = A w, and

```math
p(A,x)=\beta(x)\frac{\sqrt A-\sqrt{A_0(x)}}{A_0(x)},
\qquad \beta=\frac{Eh\sqrt\pi}{\rho(1-\xi^2)}.
```

The first-order model is

```math
\partial_t A+\partial_x Q=0,\qquad
\partial_t Q+\partial_x\left(\frac{Q^2}{A}+Ap\right)
=p\partial_x A+\Gamma_1\frac QA,\qquad \Gamma_1=2\pi Rk.
```

The second-order model adds axial diffusion and changes wall friction:

```math
\partial_t A+\partial_x Q=0,\qquad
\partial_t Q+\partial_x\left(\frac{Q^2}{A}+Ap\right)
-\partial_x\left(3\nu A\partial_x\left(\frac QA\right)\right)
=p\partial_x A+\Gamma_2\frac QA,\qquad
\Gamma_2=\frac{2\pi Rk}{1-Rk/(4\nu)}.
```

Trixi adds the divergence of the parabolic flux to the right-hand side.
The implemented flux is therefore positive: 3nu (Q_x - Q/A A_x).
With the chosen friction law, both friction and diffusion vanish continuously
at nu = 0; this case is handled without dividing by zero.

### Energy and characteristics

Define

```math
\widetilde p=\frac{\beta}{3A_0}(A^{3/2}-A_0^{3/2}),\qquad
\mathcal E_1=\frac{Q^2}{2A}+Ap-\widetilde p.
```

For smooth solutions of the second-order model,

```math
\partial_t\mathcal E_1+
\partial_x\left((\mathcal E_1+\widetilde p)w\right)
=w\partial_x(3\nu A\partial_x w)+\Gamma_2w^2.
```

With vanishing boundary energy flux and diffusion work, integration gives

```math
\frac{d}{dt}\int_0^L\mathcal E_1\,dx
=-\int_0^L3\nu A(\partial_x w)^2\,dx
+\int_0^L\Gamma_2w^2\,dx\leq0.
```

For the first-order model, omit diffusion and use Gamma1 instead of Gamma2.
At open boundaries, retain the boundary energy flux. Diffusion work need not
be non-positive pointwise.

The hyperbolic speeds are w +/- c, where c = sqrt(A P_A/rho).
For locally fixed material fields, W+ = w + 4c and W- = w - 4c are Riemann
invariants of the source-free hyperbolic part. Outflow preserves the outgoing
invariant and imposes the incoming invariant of the rest state.
It assumes subcritical axial outflow. The second-order model uses this same
hyperbolic condition together with the parabolic boundary operators.

## Two-dimensional model

The coordinates are angle theta and axial arclength s. Here A = R^2/2,
M = QRtheta = (3/4) R A wtheta, and N = Qs = A ws.
This A differs from the full cross-sectional area of the 1D model.

```math
p=P/\rho=b\frac{R-R_0}{R_0^2},\qquad
b=\frac{Eh}{\rho(1-\xi^2)},\qquad R=\sqrt{2A}.
```

Writing C(s) for curvature, the equations are

```math
\begin{aligned}
\partial_t A+\partial_\theta(M/A)+\partial_sN&=0,\\
\partial_t M+\partial_\theta\left(\frac{M^2}{2A^2}+Ap\right)
+\partial_s(MN/A)
&=p\partial_\theta A+\frac{2R}{3}C\sin\theta\frac{N^2}{A}
+2Rk\frac MA,\\
\partial_t N+\partial_\theta(MN/A^2)
+\partial_s\left(\frac{N^2}{A}-\frac{M^2}{2A^2}+Ap\right)
&=p\partial_s A-\frac{2R}{3}C\sin\theta\frac{MN}{A^2}
+Rk\frac NA.
\end{aligned}
```

### Energy and characteristics

With beta2 = b/sqrt(2), define

```math
\widetilde p=\frac{\beta_2}{3A_0}(A^{3/2}-A_0^{3/2}),\qquad
\mathcal E_2=\frac{M^2}{2A^2}+\frac{N^2}{2A}+Ap-\widetilde p
=A\left(\frac{9}{16}w_\theta^2+\frac12w_s^2+p\right)-\widetilde p.
```

Let G = E2 + ptilde - M^2/(2A^2). The energy identity is

```math
\partial_t\mathcal E_2+\partial_\theta\left(\frac{M}{A^2}G\right)
+\partial_s\left(\frac NA G\right)
=\frac94Rk w_\theta^2+Rk w_s^2\leq0.
```

Curvature sources cancel in this identity. Global energy decay additionally
requires zero or appropriately controlled boundary energy flux.

The coordinate-direction characteristic speeds are

```math
\lambda_\theta\in\left\{-\sqrt{p_A},\sqrt{p_A},M/A^2\right\},\qquad
\lambda_s\in\left\{w_s-c,w_s,w_s+c\right\}.
```

For positive compliance, A p_A > (9/8) wtheta^2 gives a positive definite
energy Hessian and strict hyperbolicity in every spatial direction.
Checking distinct eigenvalues only along the coordinate axes is insufficient
for this stronger statement; see the audit.

Axial outflow uses ws +/- 4c and preserves M/A, which remains constant across
axial acoustic waves. Characteristic reconstruction applies to boundaries
aligned with s; other normals use extrapolation.
