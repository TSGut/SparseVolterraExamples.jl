#############################
##  Numerical experiments associated with section 5.3.2
#############################
## This file contains the numerical experiments used to compare our method with
## Chebfun's Volterra integral equation implementation when high polynomial orders are required.
## We only include our own side of the implementation here. See the references in the paper for Chebfun.
##

using ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, SpecialFunctions, Plots, SparseArrays
using SparseVolterraExamples

#####
## First, we can check that the stated arctan function indeed approximates a step-like function.
## This is in Figure 5(a).
v(x,K) = atan(K*x)
    plot(x->v(x,10),0,1, legend=:true , xlabel = "x", label="k=10" ,ylabel = "u_2(x,K)",  legendfontsize=12, tickfontsize=10, thickness_scaling = 1.2 , grid=:none)
    plot!(x->v(x,50),0,1, legend=:true , xlabel = "x", label="k=50" ,ylabel = "u_2(x,K)",  legendfontsize=12, tickfontsize=10, thickness_scaling = 1.2 , grid=:none)
    plot!(x->v(x,100),0,1, legend=:true , xlabel = "x", label="k=100" , ylabel = "u_2(x,k)", legendfontsize=12, tickfontsize=10, thickness_scaling = 1.2 , grid=:none)
    plot!(x->v(x,200),0,1, legend=:bottomright , xlabel = "x", label="k=200" , ylabel = "u_2(x,k)", legendfontsize=12, tickfontsize=10, thickness_scaling = 1.2 , grid=:none)

#####
## Now we work towards Table 2 and Figure 4.
## First we define the Kernel and g(x,k) as in section 5.3.2.
gf(x,k) = k/(k^2*x^2+1)-(exp(x^2)*atan(k*x))/(2*k^2)+(exp(x^2)*x)/(2*k)-1/2*exp(x^2)*x^2*atan(k*x)
Kfun(x,y) = y*exp(x^2)

#####
## Solver function with step-by-step explanation for given k and with polynomial degree 'n'
function solveSec532(k,n,gf,Kfun)
    gF = Fun(x->gf(x,k),Jacobi(1,2, 0..1),n)             # Approximate g in the appropriate basis
    V = triVolterraFullKernelOpP01(Kfun,n,true,155)    # The following steps generate the appropriate Volterra operator
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V
        V = Conversion(Jacobi(0,0,0..1),Jacobi(1,2,0..1))[1:n,1:n]*V[1:n,1:n]
        V = Derivative(Jacobi(0,1,0..1),1)[1:n,1:n]-V
    coeff = [DirectEvalLHSP10at0(n);V[1:n-1,1:n]] \ [0;pad(gF.coefficients,n-1)] # Append evaluation and initial conditions and then solve the equation
    return Fun(Jacobi(0,1, 0..1),coeff) # returns solution Fun
end

#####
## Now we can compute and plot a specific example for k=100, plotting both the numerical approximation and the analytic solution.
## Computing errors for various k is a bit more work since automatic convergence testing is not currently implemented but a straightforward for loop can do it from here.
u = solveSec532(100,50,gf,Kfun)
plot(u, label="sparse method")
    plot!(x->v(x,100),0,1, legend=:bottomright , xlabel = "x", label="analytic" ,ylabel = "u_2(x,100)",  legendfontsize=12, tickfontsize=10, thickness_scaling = 1.2 , grid=:none)

#####
## We can also use BenchmarkTools to get a somewhat robust time estimation.
## Note however that this obviously strongly depends on the hardware you are using.
using BenchmarkTools
@benchmark solveSec532(100,50,gf,Kfun)
@benchmark solveSec532(200,300,gf,Kfun)

#####
## Now we plot the operator bandedness for the k=100 example, this is basically Figure 5(b).
function OperatorSec532(k,n,Kfun)
    V = triVolterraFullKernelOpP01(Kfun,n,true,155)    # The following steps generate the appropriate Volterra operator
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V
        V = Conversion(Jacobi(0,0,0..1),Jacobi(1,2,0..1))[1:n,1:n]*V[1:n,1:n]
        V = Derivative(Jacobi(0,1,0..1),1)[1:n,1:n]-V
        V = [DirectEvalLHSP10at0(n);V[1:n-1,1:n]] # Append evaluation and initial conditions
    return V # returns operator
end
#####
## Plots.jl's spy plot is not currently compatible with these types, so we instead convert to generic sparse for visualization.
## The exact bandedness properties depend on the chosen degrees and parameters but here is a standard example
spy(sparse(OperatorSec532(100,300,Kfun)),markersize=2.8,marker=:rect)
