#############################
##  Numerical experiments associated with section 5.4.1
#############################
## Nonlinear equations are treated by first defining an objective function, then using standard iterative methods
## such as Newton's method via NLsolve.jl. Using more sophisticated methods, linesearches etc. can sometimes result in significant speed boosts.
##

using NLsolve, ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, Plots
using SparseVolterraExamples

#####################################################
## Problem in Equation (25), u_1(x) = exp(x)
#####################################################
####
## Define g(x), the nonlinearity and the kernel
g(x) = exp(x)+1/3*x*(1-exp(3*x))
    gFun = Fun(g,Jacobi(0,1,0..1))
    f(u) = u^3
    K(x,y) = x
    n = 15 # Desired polynomial order of approximation for accuracy
####
## Define the objective function and compute the solution using a standard Newton method
objective(x) = triNonLinearVolterraObjective(x,f,K,gFun,n)
sol = nlsolve(objective,zeros(n),method=:newton,ftol=1e-8)
####
## Plot function aganinst analytic solution
plot(Fun(Jacobi(0,1,0..1),sol.zero),grid=false,legend=:topleft,xlabel="x",ylabel="u(x)",label="sparse method")
plot!(exp,0,1,grid=false,legend=:topleft,xlabel="x",ylabel="u(x)",label="analytic exp(x)")
####
## Or plot error. A simple loop over polynomial degree can then produce Figure 7(a).
plot(Fun(Jacobi(0,1,0..1),sol.zero)-Fun(x->exp(x),Jacobi(0,1,0..1)),grid=false,legend=:bottomleft,xlabel="x",ylabel="error",label=false)

#####################################################
## Problem in Equation (26), u_2(x) = sin(x)
#####################################################
####
## Define g(x), the nonlinearity and the kernel
g(x) = sin(x)+1/4*sin(x)^2-1/4*x^2
    gFun = Fun(g,Jacobi(0,1,0..1))
    f(u) = u^2
    K(x,y) = x-y
    n = 14 # Desired polynomial order of approximation for accuracy
####
## Define the objective function and compute the solution using a standard Newton method
objective(x) = triNonLinearVolterraObjective(x,f,K,gFun,n)
sol = nlsolve(objective,zeros(n),method=:newton,ftol=1e-15)
####
## Plot function aganinst analytic solution
plot(Fun(Jacobi(0,1,0..1),sol.zero),grid=false,legend=:topleft,xlabel="x",ylabel="u(x)",label="sparse method")
plot!(sin,0,1,grid=false,legend=:topleft,xlabel="x",ylabel="u(x)",label="analytic sin(x)")
####
## Or plot error. A simple loop over polynomial degree can then produce Figure 7(b).
plot(Fun(Jacobi(0,1,0..1),sol.zero)-Fun(x->sin(x),Jacobi(0,1,0..1)),grid=false,legend=:bottomleft,xlabel="x",ylabel="error",label=false)
