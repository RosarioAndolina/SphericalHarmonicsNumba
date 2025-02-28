# Spherical Harmonics Module

## Introduction

Spherical harmonics are special functions defined on the surface of a sphere. They are widely used in physics, chemistry, and engineering, particularly in solving problems with spherical symmetry, such as the Schrödinger equation for the hydrogen atom.

This Python module provides an implementation of spherical harmonics using associated Legendre polynomials. The implementation is optimized using Numba for performance and includes a comparison with SciPy's `sph_harm` function for validation.

## Theory

### Spherical Harmonics

Spherical harmonics \( Y_{l}^{m}(\theta, \phi) \) are solutions to the Laplace equation in spherical coordinates. They are defined as:

\[
Y_{l}^{m}(\theta, \phi) = \sqrt{\frac{(2l + 1)}{4\pi} \frac{(l - m)!}{(l + m)!}} P_{l}^{m}(\cos \theta) e^{im\phi}
\]

where:
- \( l \) is the degree (azimuthal quantum number),
- \( m \) is the order (magnetic quantum number),
- \( P_{l}^{m} \) is the associated Legendre polynomial,
- \( \theta \) is the polar angle,
- \( \phi \) is the azimuthal angle.

### Associated Legendre Polynomials

The associated Legendre polynomials \( P_{l}^{m}(x) \) are solutions to the associated Legendre differential equation. They are defined as:

\[
P_{l}^{m}(x) = (-1)^m (1 - x^2)^{m/2} \frac{d^m}{dx^m} P_l(x)
\]

where \( P_l(x) \) is the Legendre polynomial of degree \( l \).

The recursive algorithm used in this module computes \( P_{l}^{m}(x) \) efficiently.

## Code Overview

The module consists of the following functions:
1. `factorial(n)`: Computes the factorial of an integer.
2. `associated_legendre(l, m, x)`: Computes the associated Legendre polynomial \( P_{l}^{m}(x) \).
3. `sph_harm(m, l, phi, theta)`: Computes the spherical harmonic \( Y_{l}^{m}(\theta, \phi) \).
4. `test_specific_case(l, m, theta, phi)`: Compares the custom implementation with SciPy's `sph_harm`.

## Installation

To install the module, run the following command in the root directory of the project:

```bash
pip install .

### Usage example

```python
import numpy as np
from sph_harm_module import sph_harm

l = 2
m = 1
theta = np.pi / 4
phi = np.pi / 2

Y = sph_harm(m, l, phi, theta)
print(Y)

