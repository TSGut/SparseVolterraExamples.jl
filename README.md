[![DOI](https://zenodo.org/badge/319677914.svg)](https://zenodo.org/badge/latestdoi/319677914)

# SparseVolterraExamples.jl
This Julia package contains examples which implement the method described in https://arxiv.org/abs/2005.06081 for solving nonlinear and integro-differential Volterra equations. It additionally also contains some linear examples.

In its current state, this package serves two purposes:
- Complimentary role to the paper, providing code samples to those curious about the actual implementation.
- Reproduce the figures found in the paper by providing a documented (via comments) step-by-step walkthrough of the numerical experiments in the paper.

While the package is not meant for nor setup for user-friendly general use, the interested reader should nevertheless find the code to be straightforward to generalize.

# Installation

The following installation guide assumes a clean install of Julia Version 1.3.1 (2019-12-30), which can be found [here](https://julialang.org/downloads/oldreleases/). This is the version in which the package was written and tested.
Newer Julia releases may work if all the appropriate dependencies can be resolved.

As an unregistered Julia package, you can install this package via ```] add https://github.com/TSGut/SparseVolterraExamples.jl``` or alternatively ```Pkg.add(PackageSpec(url="https://github.com/TSGut/SparseVolterraExamples.jl"))```. Test whether everything works as intended by executing ```Pkg.test("SparseVolterraExamples")``` (this may take a few minutes depending on your hardware). Assuming this executes without errors, you can then explore the files in the examples folder.

# Using the package

The examples folder contains scripts with detailed comments. Evaluating line-by-line from the top walks the user through the process of obtaining the solutions described in the paper. Along with the comments in the src folder this should be sufficient to generalize it to other use cases. Keep in mind that to run the examples, you will have to install the packages specified via ```using``` at the top of each example and are part of the dependencies of this package.

# References

- Gutleb, T. S. (2021) ‘A Fast Sparse Spectral Method for Nonlinear Integro-Differential Volterra Equations with General Kernels’. Advances in Computational Mathematics, vol. 47, no. 3, https://doi.org/10.1007/s10444-021-09866-7.

- Gutleb, T. S., & Olver, S. (2020). 'A Sparse Spectral Method for Volterra Integral Equations Using Orthogonal Polynomials on the Triangle'. SIAM Journal on Numerical Analysis, 58(3), 1993-2018. https://doi.org/10.1137/19M1267441
