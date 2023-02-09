[![DOI](https://zenodo.org/badge/319677914.svg)](https://zenodo.org/badge/latestdoi/319677914)

# SparseVolterraExamples.jl
This Julia package contains examples which implement the method described in https://arxiv.org/abs/2005.06081 for solving nonlinear and integro-differential Volterra equations. It additionally also contains some linear examples.

In its current state, this package serves two purposes:
- Complimentary role to the paper, providing code samples to those curious about the actual implementation.
- Reproduce the figures found in the paper by providing a documented (via comments) step-by-step walkthrough of the numerical experiments in the paper.

While the package is not meant for nor setup for user-friendly general use, the interested reader should nevertheless find the code to be straightforward to generalize.

# Installation

The current version has been tested with Julia 1.8.5. Newer Julia releases may work if all the appropriate dependencies can be resolved.
For the version used in the papers which includes tests and has been tested up to Julia v.1.5.3, see [v0.1.1](https://github.com/TSGut/SparseVolterraExamples.jl/releases/tag/v0.1.1)

As an unregistered Julia package, you can install this package via ```] add https://github.com/TSGut/SparseVolterraExamples.jl``` or alternatively ```Pkg.add(PackageSpec(url="https://github.com/TSGut/SparseVolterraExamples.jl"))```. Then explore the files in the examples folder.

# Using the package

The examples folder contains scripts with detailed comments. Evaluating line-by-line from the top walks the user through the process of obtaining the solutions described in the paper. Along with the comments in the src folder this should be sufficient to generalize it to other use cases. Keep in mind that to run the examples, you will have to install the packages specified via ```using``` at the top of each example and are part of the dependencies of this package.

# References

- Gutleb, T. S. (2021) ‘A Fast Sparse Spectral Method for Nonlinear Integro-Differential Volterra Equations with General Kernels’. Advances in Computational Mathematics, vol. 47, no. 3, https://doi.org/10.1007/s10444-021-09866-7.

- Gutleb, T. S., & Olver, S. (2020). 'A Sparse Spectral Method for Volterra Integral Equations Using Orthogonal Polynomials on the Triangle'. SIAM Journal on Numerical Analysis, 58(3), 1993-2018. https://doi.org/10.1137/19M1267441
