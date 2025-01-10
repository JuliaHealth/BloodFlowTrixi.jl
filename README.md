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

- **2D Blood Flow Model**: The 2D model extends the Navier-Stokes equations under the thin-artery assumption, allowing for simulations in complex arterial geometries using curvilinear coordinates. It captures both longitudinal and angular dynamics, making it more accurate than classical 1D models while being less computationally expensive than full 3D models.

  This model is described in detail in:  
  **[[Article 2D](https://hal.science/hal-04700161v1)]**

Both models were designed to be used with **[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)**, a flexible and high-performance framework for solving systems of conservation laws using the Discontinuous Galerkin (DG) method.

## Features

- **1D and 2D models** for arterial blood flow.
- Derived from the Navier-Stokes equations with appropriate assumptions for compliant arteries.
- To be used with **Trixi.jl** for DG-based numerical simulations.

## 
- Support for curvilinear geometries and compliant wall dynamics.
- Validated through peer-reviewed research articles.

## Installation

To install **BloodFlowTrixi.jl**, use the following commands in Julia:

```bash
julia> ]
pkg> add Trixi
pkg> add https://github.com/your-repo/BloodFlowTrixi.jl
```

## License

This package is licensed under the MIT license.

## Acknowledgments

This package was developed as part of my PhD research in applied mathematics, focusing on mathematical modeling and numerical simulation of blood flow in arteries. Special thanks to the developers of **Trixi.jl**, whose framework was invaluable in implementing and testing these models.

