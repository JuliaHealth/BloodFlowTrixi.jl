# BloodFlowTrixi.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yolhan83.github.io/BloodFlowTrixi.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yolhan83.github.io/BloodFlowTrixi.jl/dev/)
[![Build Status](https://github.com/yolhan83/BloodFlowTrixi.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yolhan83/BloodFlowTrixi.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**BloodFlowTrixi.jl** is a Julia package that implements one-dimensional (1D) and two-dimensional (2D) blood flow models for arterial circulation. These models are derived from the Navier-Stokes equations and were developed as part of my PhD research in applied mathematics, focusing on cardiovascular pathologies such as aneurysms and stenoses.

## Description

This package provides:

- **1D Blood Flow Model**: This model describes blood flow along compliant arteries in a single spatial dimension. It was derived under the assumption of axisymmetric flow and accounts for arterial compliance, inertia, and frictional losses.
  
  More details about this model can be found in my corresponding publication:  
  **[[Article 1D](https://hal.science/hal-04676130v1)]**
```math
\left\{\begin{aligned}
  \frac{\partial a}{\partial t} + \frac{\partial}{\partial x}(Q) &= 0 \\
  \frac{\partial Q}{\partial t} + \frac{\partial}{\partial x}\left(\frac{Q^2}{A} + A P(a)\right) &= P(a) \frac{\partial A}{\partial x} - 2 \pi R k \frac Q {A}\\
  P(a) &= P_{ext} + \frac{Eh\sqrt{\pi}}{1-\xi^2}\frac{\sqrt{A} - \sqrt{A_0}}{A_0} \\
  R &= \sqrt{\frac{A}{\pi}}
\end{aligned}\right.
```

- **2D Blood Flow Model**: The 2D model extends the Navier-Stokes equations under the thin-artery assumption, allowing for simulations in complex arterial geometries using curvilinear coordinates. It captures both longitudinal and angular dynamics, making it more accurate than classical 1D models while being less computationally expensive than full 3D models.

  This model is described in detail in:  
  **[[Article 2D](https://hal.science/hal-04700161v1)]**
```math
\left\{\begin{aligned}
    \frac{\partial a}{\partial t} + \frac{\partial}{\partial \theta}\left( \frac{Q_{R\theta}}{A} \right) + \frac{\partial}{\partial s}(Q_s) &= 0 \\
    \frac{\partial Q_{R\theta}}{\partial t} + \frac{\partial}{\partial \theta}\left(\frac{Q_{R\theta}^2}{2A^2} + A P(a)\right) + \frac{\partial}{\partial s}\left( \frac{Q_{R\theta}Q_s}{A} \right) &= P(a) \frac{\partial A}{\partial \theta} - 2 R k \frac{Q_{R\theta}}{A} + \frac{2R}{3} \mathcal{C}\sin \theta \frac{Q_s^2}{A} \\
    \frac{\partial Q_{s}}{\partial t} + \frac{\partial}{\partial \theta}\left(\frac{Q_{R\theta} Q_s}{A^2} \right) + \frac{\partial}{\partial s}\left( \frac{Q_s^2}{A} - \frac{Q_{R\theta}^2}{2A^2} + A P(a) \right) &= P(a) \frac{\partial A}{\partial s} - R k \frac{Q_s}{A} - \frac{2R}{3} \mathcal{C}\sin \theta \frac{Q_s Q_{R\theta}}{A^2} \\
    P(a) &= P_{ext} + \frac{Eh}{\sqrt{2}\left(1-\xi^2\right)}\frac{\sqrt{A} - \sqrt{A_0}}{A_0} \\
    R &= \sqrt{2A}
\end{aligned}\right.
```

Both models were designed to be used with **[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)**, a flexible and high-performance framework for solving systems of conservation laws using the Discontinuous Galerkin (DG) method.

## Features

- **1D and 2D models** for arterial blood flow.
- Derived from the Navier-Stokes equations with appropriate assumptions for compliant arteries.
- To be used with **Trixi.jl** for DG-based numerical simulations.
- Support for curvilinear geometries and compliant wall dynamics.

## Installation

To install **BloodFlowTrixi.jl**, use the following commands in Julia:

```bash
julia> ]
pkg> add Trixi
pkg> add BloodFlowTrixi
```

## Future Plans

**short term**
- Add second order 1D model.
- Add proper tests for 1D and 2D models.
- Add 3D representations of the solutions for 1D and 2D models.
- Design easy to use interfaces for users to define their own initial and boundary conditions and source terms.

**long term**
- Add 3D fluid-structure interaction models for complex arterial geometries.
- Design support for artery networks and simulate vascular networks using the 2D and 1D model.
- Autodiff support for 1D and 2D models for parameter optimization.

## License

This package is licensed under the MIT license.

## Acknowledgments

This package was developed as part of my PhD research in applied mathematics, focusing on mathematical modeling and numerical simulation of blood flow in arteries. Special thanks to the developers of **Trixi.jl**, whose framework was invaluable in implementing and testing these models.

