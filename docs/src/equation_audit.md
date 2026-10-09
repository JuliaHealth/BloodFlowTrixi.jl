# Equation verification

This audit compares the reduced PDEs, pressure laws, source coefficients,
characteristic spectra and energy identities with the published papers. It does
not constitute an independent proof of the complete asymptotic derivation from
three-dimensional Navier-Stokes equations.

## References and equations used

- [1D paper, DOI 10.4236/jamp.2025.1310198](https://content.scirp.org/pdf/jamp_1724375.pdf):
  first-order model (31), printed p. 3489; viscous models (5) and (38), pp. 3481 and 3494.
- [2D paper, DOI 10.4236/jamp.2025.1311220](https://content.scirp.org/pdf/jamp_1724376.pdf):
  model (5)/(41), pp. 3932 and 3945; energy and directional eigenvalues, pp. 3946-3949.

The main reduced equations agree with their pressure laws, characteristic spectra
and local energy identities when using the corrections below. The code uses
physical transmural pressure P; the derivation variable p is P/rho.
The canonical equations are given in [Mathematics](math.md).

## Corrections needed in the 1D paper

1. **Theorem 1, item 5, p. 3489 and its proof p. 3491:** the global energy
   formula contains the second-order denominator `1 - epsilon Rk/(4nu0)`,
   although model (31) and local identity (34) use first-order friction.
   With zero boundary flux, the correct identity is

   ```math
   \frac{d}{dt}\int\mathcal E_1\,dx=\int 2\pi Rk\,w^2\,dx.
   ```

2. **Theorem 2, items 1 and 3, p. 3494, and its proof p. 3495:** the local
   head and energy equations omit the denominator present in model (38).
   Write `Gamma2 = 2pi Rk/(1 - epsilon Rk/(4nu0))`; their friction terms must
   then be `Gamma2 * w/A` and `Gamma2 * w^2`, respectively. Introductory
   equation (5) includes this denominator correctly in physical variables.

3. **Global energy formulas:** if Rk varies with position, its coefficient
   belongs inside the integral. It cannot multiply a squared L2 norm as a
   single constant. The code's choice k = -11nu/R makes Rk constant.

These are inconsistencies in the stated energy formulas, rather than changes
to the principal models (5), (31) and (38).

## Clarification needed in the 2D paper

Theorem 1, item 1, p. 3946 checks distinct eigenvalues for the coordinate
matrices, using `A p_A - (9/8) wtheta^2 != 0`. This does not suffice for strict
hyperbolicity in **every** spatial direction.

For example, take A = 1, M = 2, N = 0, p_A = 1. The stated condition is
1 - 4 = -3 != 0, and both coordinate matrices have distinct spectra.
For the unnormalized normal n = (1, sqrt(3)), their linear combination is

```math
H_n=\begin{pmatrix}
-2&1&\sqrt3\\
-3&2&2\sqrt3\\
5\sqrt3&-2\sqrt3&2
\end{pmatrix}.
```

Its characteristic polynomial is (lambda - 2)^2 (lambda + 2), and
rank(Hn - 2I) = 2: the repeated eigenvalue has only one independent
eigenvector. The system is therefore not strictly hyperbolic in this direction.
The sufficient condition `A p_A > (9/8) wtheta^2` also makes the energy Hessian
positive definite. If the paper means only coordinatewise strict hyperbolicity,
that restriction should be stated explicitly.

The energy coefficient is **9/16** in front of A wtheta^2: the 9/8 in the
numerator defining the kinetic head on p. 3946 is itself divided by two.
The principal equations use angular friction 2Rk M/A and axial curvature
-(2R/3) C sin(theta) MN/A^2.

Global energy decay on p. 3949 additionally requires the boundary energy flux
to vanish or have an appropriate sign. Local dissipation alone does not imply
global decay for arbitrary inflow/outflow data.

## Geometry of the constant-curvature example

The coordinate Jacobian in the 2D derivation contains
`beta = 1 - r C(s) cos(theta)` (p. 3934). A regular tubular parameterization
requires R C < 1 throughout the vessel; the thin-vessel approximation additionally
uses small geometric corrections. The simple example has C = 1 and R0 = 2,
so it exceeds this geometric condition. It is retained as requested as a numerical
stress test of the reduced PDE, rather than a validation within the geometric
assumptions of the derivation. Numerical stability does not establish geometric
admissibility.

## Implementation and regression checks

The implementation now uses P/rho in conservative and nonconservative fluxes,
energy and wave speeds, preserving physical pressure in the public API and
prescribed-pressure boundaries. Inverse pressure functions keep their physical
pressure convention; characteristic reconstruction includes density.
Angular friction uses the published factor 2, and axial curvature uses A^2.
The zero-viscosity limit of second-order friction is finite.

Tests differentiate the implemented fluxes and nonconservative jumps to recover
the continuous operators. They compare their eigenvalues with the paper spectra,
check the differential energy-flux identity including material-field columns,
and compare source work with analytic friction dissipation for densities 1 and
1.7 and three area ratios. Characteristic boundaries are also exercised with
both densities in Float32 and Float64. A parabolic RHS test checks energy
dissipation; manufactured solutions check forcing against the full PDE.

The two convergence forcings were corrected: the 1D momentum residual lacked a
factor 2, and the 2D forcing used the wrong reference-area index and returned
four entries for a five-variable system. These prescribed forcing terms replace
physical curvature/friction, rather than adding to the physical sources.

The 2D visualization conversion now retains `wtheta = 4M/(3RA)` instead of
overwriting it with M/A, and 1D visualization returns its documented physical
pressure field.
