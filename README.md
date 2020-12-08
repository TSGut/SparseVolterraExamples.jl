# SparseVolterraExamples.jl
This Julia package containts examples which implement the method described in https://arxiv.org/abs/2005.06081 for solving nonlinear and integro-differential Volterra equations. It additionally also contains some linear examples.

In its current state, this package serves two purposes: 
- Complimentary role to the paper, providing code samples to those curious about the actual implementation.
- Reproduce the figures found in the paper by providing a documented (via comments) reproduction of each of the numerical experiments in the paper.

While the package is not meant for, nor setup for user-friendly standalone use, the interested reader should nevertheless find the code to be mostly straightforward to generalize to their use-case.

# References

- Gutleb, T. S. (2020). A fast sparse spectral method for nonlinear integro-differential Volterra equations with general kernels. arXiv:2005.06081.

- Gutleb, T. S., & Olver, S. (2020). A sparse spectral method for Volterra integral equations using orthogonal polynomials on the triangle. SIAM Journal on Numerical Analysis, 58(3), 1993-2018.
