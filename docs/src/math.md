# Summary of the 1D and 2D Blood Flow Models and Numerical Approximation

## 1D Blood Flow Model

### Model Assumptions (1D)

1. **H1**: Thin wall and plane stresses — the artery wall thickness is constant and small enough to allow a shell-type representation.
2. **H2**: Radial displacement — the artery displacements occur only in the radial direction.
3. **H3**: Small deformation gradients and linear elastic behavior — the artery behaves as a linear elastic solid.
4. **H4**: Incompressibility — the wall tissue is incompressible.
5. **H5**: Dominance of circumferential stresses — axial stresses are negligible compared to circumferential ones.

### Asymptotic Analysis (1D)

The model is derived from the incompressible Navier-Stokes equations in cylindrical coordinates under the **thin-artery assumption** (small ratio of radius to length). A section-averaged approach is applied to reduce the three-dimensional problem to a one-dimensional model.

#### Non-Dimensionalization

Dimensionless variables are introduced using the following scaling:

- Time: \( \tilde{t} = \frac{t}{T} \)
- Space: \( \tilde{x}, \tilde{r} \)
- Pressure: \( \tilde{p} \)
- Velocity: \( \tilde{u}_x, \tilde{u}_r \)

A small parameter \( \epsilon \) representing the ratio of radius to length is used for asymptotic expansion.

### 1D Model Equations

The derived one-dimensional section-averaged model is given by:

```math
\[
\begin{aligned}
& \partial_t A + \partial_x Q = 0, \\
& \partial_t Q + \partial_x \left( \frac{Q^2}{A} + \frac{1}{\rho} P(A, x) \right) - \partial_x \left( 3 \nu A \partial_x \left( \frac{Q}{A} \right) \right) = \frac{1}{\rho} P(A, x) \partial_x A + \frac{2 \pi R k}{1 - R k / (4 \nu)} \frac{Q}{A},
\end{aligned}
\]
```

where \( A \) is the cross-sectional area, \( Q \) is the flow rate, and \( P(A, x) \) is the pressure given by:

```math
\[
P(A, x) = b(x) \frac{\sqrt{A} - \sqrt{A_0}}{A_0},
\]
```

with \( b(x) \) representing the artery's elastic properties.

### Energy Consistency Theorem (1D)

Let \((A, u_x)\) be a smooth solution of the system, where \( u_x = \frac{Q}{A} \) is the mean velocity. The total energy \( E \) is defined by:

```math
\[
E = A \left( \frac{u_x^2}{2} + \frac{1}{\rho} P(A, x) \right) - \tilde{P},
\]
```

where \( \tilde{P} \) is such that \( \tilde{P}'(A) = A P'(A) \).

The entropy pair \((E, E + \tilde{P})\) satisfies the following entropy relation:

```math
\[
\partial_t E + \partial_x \left( \left( E + \tilde{P} \right) u_x \right) = 3 \nu A (\partial_x u_x)^2 + \frac{2 \pi R k}{1 - R k / (4 \nu)} u_x^2 \leq 0.
\]
```

This shows that the total energy decreases over time, accounting for viscous dissipation and friction effects.

---

## 2D Blood Flow Model

### Model Assumptions (2D)

1. **H1**: Thin wall and plane stresses — the vessel wall thickness is assumed constant and sufficiently thin to allow a shell-type representation.
2. **H2**: Radial displacement — displacements of the artery occur only in the radial direction.
3. **H3**: Small deformation gradients and linear elastic behavior — the artery wall behaves as a linear elastic solid.
4. **H4**: Incompressibility — the artery wall tissue is incompressible.
5. **H5**: Dominance of circumferential stresses — stresses acting in the axial direction can be neglected compared to those in the circumferential direction.

### Asymptotic Analysis (2D)

Under the **thin-artery assumption** (small ratio between the radius and length), an asymptotic expansion of the Navier-Stokes equations was performed to first order in terms of a parameter \( \epsilon \). This leads to a radially averaged two-dimensional model while neglecting higher-order terms in \( \epsilon \).

The dimensional equations are rewritten in non-dimensional form using appropriately scaled variables:

- Time: \( \tilde{t} = \frac{t}{T} \)
- Space: \( \tilde{s}, \tilde{r}, \tilde{\theta} \)
- Pressure: \( \tilde{p} \)
- Velocity: \( \tilde{u}_s, \tilde{u}_r, \tilde{u}_\theta \)

### 2D Model Equations

The radially averaged two-dimensional system is given by:

```math
\[
\begin{aligned}
& \partial_t A + \partial_\theta \left( \frac{Q_{R\theta}}{A} \right) + \partial_s(Q_s) = 0, \\
& \partial_t (Q_{R\theta}) + \partial_\theta \left( \frac{Q_{R\theta}^2}{2A^2} + Ap \right) + \partial_s \left( \frac{Q_{R\theta}Q_s}{A} \right) = \frac{2R}{3} C \sin \theta \frac{Q_s^2}{A} + \frac{2Rk Q_{R\theta}}{A} + \partial_\theta Ap, \\
& \partial_t (Q_s) + \partial_\theta \left( \frac{Q_s Q_{R\theta}}{A^2} \right) + \partial_s \left( \frac{Q_s^2}{A} - \frac{Q_{R\theta}^2}{2A^2} + Ap \right) = - \frac{2R}{3} C \sin \theta \frac{Q_s Q_{R\theta}}{A^2} + \frac{kR Q_s}{A} + \partial_s Ap, \\
& p = p_{\text{ext}} + b \frac{R - R_0}{R_0^2}.
\end{aligned}
\]
```

### Energy Consistency Theorem (2D)

Let \((A, u_\theta, u_s)\) be a smooth solution of the 2D system. The total energy \( E \) is defined as:

```math
\[
E = A \left( \psi + p \right) - \tilde{p},
\]
```

where \( \psi = \frac{9}{16} u_\theta^2 + \frac{1}{2} u_s^2 \) is the total head, and \( \tilde{p} \) is such that \( \tilde{p}'(A) = A p'(A) \).

The entropy pair \((E, E + \tilde{p})\) satisfies the following entropy relation:

```math
\[
\partial_t E + \nabla_{\theta, s} \cdot \left( \begin{bmatrix}
\frac{3}{2} \frac{u_\theta}{R} \\
u_s
\end{bmatrix} (E + \tilde{p} - \frac{9}{16} A u_\theta^2) \right) = \frac{9}{4} R k u_\theta^2 + k R u_s^2.
\]
```

This relation shows that the total energy of the system decreases over time, accounting for friction effects (terms in \( k \)).

