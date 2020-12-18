#############################
##  Numerical experiments associated with section 5.5.2
#############################
## In general, VIDEs are implemented as a combination of the approaches for nonlinear and integro-differential equations.
## Whether this is more efficient than standard methods once again depends on the sparsity of the resulting operator
## and the polynomial order required to resolve the solution.
##

using NLsolve, ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, Plots
using SparseVolterraExamples

#####################################################
## Problem in Equation (28), u_2(x) = tan(x)
#####################################################

####
## Some initialization
DD = 1                      # Highest derivative order of the problem
    f(u) = u^2              # This is the nonlinearity
    lincomb(x) = sin(x)     # In lincomb we store potentially appearing additional linear terms. In this case this is sin(x) as seen in the example.
    g(x) = x+cos(x)-tan(x)+tan(x)^2   # the g function from the problem statement
    K(x,y) = 1                        # the kernel is constant in this example
    conditions = [0]                  # Vector of length DD containing the conditions that uniquely solve the differential equation
    n = 28                            # Set some desired polynomial degree of approximation.
####
## Now we can set up the appropriate objective function
objective(x) = triNonLinearIntegroDiffVolterraObjective(x,f,1,g,n,DD,conditions,lincomb)

####
## Solve using basic Newton iteration with a zero vector guess
sol = nlsolve(objective,zeros(n),method=:newton,ftol=1e-15)

####
## We can compare the numerical solution against the analytic one.
plot(Fun(Jacobi(0,1,0..1),sol.zero),label="sparse method solution",grid=false,legend=:topleft,ylabel="u(x)",xlabel="x")
plot!(x->tan(x),0,1,label="analytic tan(x)",grid=false,legend=:topleft,ylabel="u(x)",xlabel="x")

####
## Or just plot the error on (0,1). A straightforward loop over the polynomial order and a maximum function over this difference then leads to Figure 8(b).
plot(Fun(Jacobi(0,1,0..1),sol.zero)-Fun(x->tan(x),0..1),label=false,grid=false,legend=:topleft,ylabel="error",xlabel="x")
