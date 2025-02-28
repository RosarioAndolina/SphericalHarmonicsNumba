import numpy as np
from numba import njit
from scipy.special import sph_harm as sph_harm_scipy

@njit
def factorial(n):
    """
    Compute the factorial of a non-negative integer using a loop.
    
    Parameters:
        n (int): The integer for which to compute the factorial.
    
    Returns:
        int: The factorial of n. Returns 0 if n is negative.
    """
    if n < 0:
        return 0
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result

@njit
def associated_legendre(l, m, x):
    """
    Compute the associated Legendre polynomial P_{l,m}(x) using a recursive algorithm.
    
    Parameters:
        l (int): Degree of the polynomial.
        m (int): Order of the polynomial.
        x (float): Point at which to evaluate the polynomial.
    
    Returns:
        float: The value of the associated Legendre polynomial P_{l,m}(x).
    """
    if m < 0:
        # Handle negative m using symmetry properties
        m = -m
        factor = (-1)**m * factorial(l - m) / factorial(l + m)
        return factor * associated_legendre(l, m, x)
    
    # Return 0 if m > l or |x| > 1 (invalid cases)
    if m > l or abs(x) > 1.0:
        return 0.0
    
    # Compute P_{m,m}(x)
    pmm = 1.0
    if m > 0:
        somx2 = np.sqrt((1.0 - x) * (1.0 + x))
        fact = 1.0
        for i in range(1, m + 1):
            pmm *= -fact * somx2
            fact += 2.0
    
    # If l == m, return P_{m,m}(x)
    if l == m:
        return pmm
    else:
        # Compute P_{m+1,m}(x)
        pmmp1 = x * (2 * m + 1) * pmm
        if l == m + 1:
            return pmmp1
        else:
            # Use recursion to compute P_{l,m}(x) for l > m + 1
            for ll in range(m + 2, l + 1):
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
                pmm = pmmp1
                pmmp1 = pll
            return pmmp1

@njit
def sph_harm(m, l, phi, theta):
    """
    Compute the spherical harmonic Y_{l,m}(theta, phi) using the associated Legendre polynomial.
    
    Parameters:
        m (int): Order of the spherical harmonic.
        l (int): Degree of the spherical harmonic.
        phi (float): Azimuthal angle (in radians).
        theta (float): Polar angle (in radians).
    
    Returns:
        complex: The value of the spherical harmonic Y_{l,m}(theta, phi).
    """
    if m < 0:
        # Handle negative m using symmetry properties
        Y = sph_harm(-m, l, phi, theta)
        return (-1)**m * np.conj(Y)
    
    # Normalization factor
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - m) / factorial(l + m))
    
    # Compute P_{l,m}(cos(theta))
    x = np.cos(theta)
    P_lm = associated_legendre(l, m, x)
    
    # Combine normalization, Legendre polynomial, and phase factor
    Y = norm * P_lm * np.exp(1j * m * phi)
    
    return Y

# Test function to compare custom implementation with SciPy's sph_harm
def test_specific_case(m, l, phi, theta):
    """
    Test the custom spherical harmonic implementation against SciPy's implementation.
    
    Parameters:
        m (int): Order of the spherical harmonic.
        l (int): Degree of the spherical harmonic.
        phi (float): Azimuthal angle (in radians).
        theta (float): Polar angle (in radians).
    """
    Y_custom = sph_harm(m, l, phi, theta)
    Y_scipy = sph_harm_scipy(m, l, phi, theta)
    
    print(f"Custom: {Y_custom}")
    print(f"SciPy:  {Y_scipy}")
    assert np.isclose(Y_custom, Y_scipy, atol=1e-10), "Error in the specific case!"


if __name__ == "__main__":
    # Test cases
    l = 2
    m = 1
    theta = 0.0
    phi = 0.0

    print(f'\nl = {l}, m = {m}, theta = {theta}, phi = {phi}')
    test_specific_case(m, l, phi, theta)

    theta = np.pi / 4
    phi = np.pi / 2

    print(f'\nl = {l}, m = {m}, theta = {theta}, phi = {phi}')
    test_specific_case(m, l, phi, theta)
