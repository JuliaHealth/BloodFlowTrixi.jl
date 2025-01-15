# 1D and 2D Mathematical Models for Blood Flow

## 1D Model
The 1D model is based on a **cross-sectional integration of the Navier-Stokes equations** under the assumption of **incompressible flow** in thin arteries. This model is particularly suitable for global studies of the arterial network, where the geometry is approximately linear or weakly curved.

### Assumptions and Simplifications
- The flow is considered **incompressible**.
- The artery is modeled as a cylindrical tube with a cross-section varying with pressure.
- A parabolic velocity profile is assumed, enabling **averaging** over the artery's cross-section.

### Main Equations
The derived equations form a system of **hyperbolic partial differential equations** describing mass and momentum conservation:

1. **Mass conservation**:
```math
∂_t A + ∂_x Q = 0
```

2. **Momentum conservation**:
```math
∂_t Q + ∂_x \left( \frac{Q^2}{A} + \frac{1}{\rho} A P(A, x) \right) - \partial_x \left( 3\nu A \partial_x\left(\frac{Q}{A}\right) \right) = \frac{1}{\rho} P(A, x) ∂_x A - \frac{2\pi R K}{1-\frac{Rk}{4\nu}} \frac{Q}{A}
```

### Energy and Entropy Relation of the 1D Model
The energy associated with the system is given by:
```math
E(t, x) = \frac{A u_x^2}{2} + \frac{1}{\rho} A P(A, x) - \frac{\beta(x)}{3 \rho A_0(x)} A^{3/2}
```

The entropy relation verified by this energy is:
```math
∂_t E + ∂_x \left( \left( E + \frac{\beta(x)}{3 \rho A_0(x)} A^{3/2} \right) u_x \right) = ∂_x \left( 3 \nu A ∂_x \left( \frac{Q}{A} \right) \right) u_x + \frac{2 \pi R k}{1 - R k / 4 \nu} u_x^2 ≤ 0
```

Under null boundary conditions:
```math
∂_t \left( \int_0^L E \, dx \right) = - 3 \nu \int_0^L A (∂_x u_x)^2 \, dx - \frac{2 \pi R k}{1 - R k / 4 \nu} \int_0^L u_x^2 \, dx < 0
```

---

## 2D Model
The 2D model is derived from a **radial integration of the Navier-Stokes equations**, enabling better representation of local effects in complex geometric configurations, such as **arterial bifurcations** and **severe aneurysms**.

### Assumptions and Simplifications
- The flow is assumed **incompressible**.
- The artery geometry is described using a curvilinear coordinate system (\( s, \theta \)).
- The velocity profile is obtained without relying on a specific ansatz.

### Main Equations
1. **Mass conservation**:
```math
∂_t A + ∂_θ \left( \frac{Q_{Rθ}}{A} \right) + ∂_s(Q_s) = 0
```

2. **Momentum conservation (radial and axial components)**:
```math
∂_t (Q_{Rθ}) + ∂_θ \left( \frac{Q_{Rθ}^2}{2 A^2} + A P \right) + ∂_s \left( \frac{Q_{Rθ} Q_s}{A} \right) = \frac{2 R}{3} C \sin θ \frac{Q_s^2}{A} + \frac{2 R k Q_{Rθ}}{A} + P∂_θ (A)
```
```math
∂_t (Q_s) + ∂_θ \left( \frac{Q_s Q_{Rθ}}{A^2} \right) + ∂_s \left( \frac{Q_s^2}{A} - \frac{Q_{Rθ}^2}{2 A^2} + A P \right) = - \frac{2 R}{3} C \sin θ \frac{Q_{Rθ} Q_s}{A^2} + \frac{k R Q_s}{A} + P∂_s (A)
```

### Energy and Entropy Relation of the 2D Model
The energy associated with the system is given by:
```math
E(t, θ, s) = A \left( \frac{9}{8} u_θ^2 + \frac{u_s^2}{2} + p \right) - \tilde{p}
```

The corresponding entropy relation is:
```math
∂_t E + ∂_θ \left( \frac{3}{2} \frac{u_θ}{R} \left( E + \tilde{p} - \frac{9}{16} A u_θ^2 \right) \right) + ∂_s \left( u_s \left( E + \tilde{p} - \frac{9}{16} A u_θ^2 \right) \right) = \frac{9}{4} R k u_θ^2 + k R u_s^2 ≤ 0
```

This relation ensures that the energy locally decreases over time, guaranteeing the **stability** of the model.

---

## Comparison of 1D and 2D Models
- **1D Model**:
  - Fast and efficient for global simulations of large arterial networks.
  - Well-suited for simple or weakly curved geometries.
  - Very low computational cost.
- **2D Model**:
  - More accurate for complex geometries (bifurcations, aneurysms).
  - Better captures local effects and fluid-structure interactions.
  - Moderate computational cost compared to three-dimensional models (3D NS-FSI).

The combined use of these two models provides an **efficient alternative to 3D simulations**, offering a good compromise between accuracy and computational cost.
