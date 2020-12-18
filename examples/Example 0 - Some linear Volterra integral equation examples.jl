#############################
##  These numerical experiments, included for the sake of completion,
##  showcase some numerical solutions to LINEAR Volterra integral equations
#############################

using ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, SpecialFunctions, Plots
using SparseVolterraExamples

#####################################################
## Solving Equation 18 in https://doi.org/10.1137/19M1267441
## This is a linear Volterra integral equation of FIRST kind.
#####################################################
####
## Set up g(x) and kernel K(x,y)
g1(x)= exp(-x)/4+1/4*exp(x)*(-1+2*x)
K(x,y) = exp(y-x)
####
## For linear problems the solver is wrapped into convenient functions. N is polynomial order of approximation.
N = 20
u = triVolterraEQ1FullKernelSolver(g1,K,N,true)
####
## Plot the computed and analytic solution
plot(u,grid=false,ylabel="u(x)",xlabel="x",label="sparse method")
plot!(x->x*exp(x),0,1,grid=false,ylabel="u(x)",xlabel="x",label="analytic")
####
## Plot error
plot(triVolterraEQ1FullKernelSolver(g1,K,N,true)-Fun(x->x*exp(x),Jacobi(0,1,0..1)),ylabel="error",xlabel="x",grid=false,label=false)

#####################################################
## The following is a simple linear second kind Volterra integral equation, not appearing in either paper.
## It uses the monomial coefficient vector input kernel implementation and limits of integration to 1-x instead and has analytically known solution sin(-2*pi*x).
#####################################################
####
## Set up g(x) and kernel K(x,y)
g2 = Fun(x->(x-1)/(2*pi)*cos(-2*pi*x)+(4*pi^2+1)/(4*pi^2)*sin(2*pi*(1-x)), Jacobi(0,1,0..1))
K2 = [0.,1.,0.,0.] # this is equivalent to K(x,y)=y
####
## Compute and plot solutions
N = 20
u = triVolterraEQ2Solver(g2,K2,N,false) # note the false indicating that we are not flipping the integration bounds, i.e. we are integrating from 0 to 1-x
plot(u,grid=false,ylabel="u(x)",xlabel="x",label="sparse method")
plot!(x->sin(-2*pi*x),0,1,grid=false,ylabel="u(x)",xlabel="x",label="analytic")
####
## Plot errors
plot(triVolterraEQ2Solver(g2,K2,N,false)-Fun(x->sin(-2*pi*x),Jacobi(0,1,0..1)),grid=false,ylabel="error",xlabel="x",label=false)
####
## The above is equivalent to the following much simpler input, which uses the full kernel input as a function
u2 = triVolterraEQ2FullKernelSolver(g2,(x,y)->y,N,false)
####
## Verify the two input methods give the same answer
plot(x->abs(u2(x)-u(x)),0,1,grid=false,ylabel="difference",xlabel="x",label=false)
