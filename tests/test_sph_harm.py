import numpy as np
from SphericalHarmonicsNumba.sph_harm import sph_harm, sph_harm_scipy

def test_sph_harm():
    """
    Test the custom spherical harmonic implementation against SciPy's implementation.
    """
    # Test case 1
    l = 2
    m = 1
    theta = 0.0
    phi = 0.0
    Y_custom = sph_harm(m, l, phi, theta)
    Y_scipy = sph_harm_scipy(m, l, phi, theta)
    assert np.isclose(Y_custom, Y_scipy, atol=1e-10), "Error in the specific case!"

    # Test case 2
    theta = np.pi / 4
    phi = np.pi / 2
    Y_custom = sph_harm(m, l, phi, theta)
    Y_scipy = sph_harm_scipy(m, l, phi, theta)
    assert np.isclose(Y_custom, Y_scipy, atol=1e-10), "Error in the specific case!"