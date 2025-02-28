from setuptools import setup, find_packages

setup(
    name="SphericalHarmonicsNumba",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "numba",
    ],
    author="Rosario Andolina",
    author_email="andolinarosario@gmail.com",
    description="Numba optimized implementation of spherical harmonics",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/RosarioAndolina/SphericalHarmonicsNumba.git",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GPL3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

