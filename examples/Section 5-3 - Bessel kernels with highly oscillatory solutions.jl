#############################
##  Numerical experiments associated with section 5.3
#############################
## This file contains the numerical experiments used for the highly oscillatory Bessel kernel example,
## which was originally featured by Nick Hale in https://doi.org/10.1093/imanum/dry042
##

using ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, SpecialFunctions, Plots
using SparseVolterraExamples

####
## First we define the Kernel and g(x) as in section 5.3.
## As mentioned in the paper, it is sensible to use pre-existing implementations to compute the Bessel functions
Kfun(x,y)=besselj(μ,ω*(x-y))
g(x) = besselj(μ+ν,ω*x)+1/(2*x^2)*((ν-1)*(ν-2)*besselj(ν-1,ω*x)+(ν+1)*(ν+2)*besselj(ν+1,ω*x))

#####
## Solver function with step-by-step explanation for given polynomial degree 'n'
## Note the factor 10^(-3) placed as discussed in the paper.
function solveSec54(n,g,Kfun,ν,μ,ω)
    gF = Fun(x->g(x),Jacobi(2,3, 0..1))
    V = triVolterraFullKernelOpP01(Kfun,n+2,true,339)
        V = reflectPabtoPba(n+2)*WLoweringP01P00(n+2)*V[1:n+2,1:n+2]
        V = Conversion(Jacobi(0,0,0..1),Jacobi(2,3,0..1))[1:n+2,1:n+2]*V[1:n+2,1:n+2]
        V = 10^(-3)*Derivative(Jacobi(0,1,0..1),2)[1:n+2,1:n+2]+ω^2*Conversion(Jacobi(0,1,0..1),Jacobi(2,3,0..1))[1:n+2,1:n+2]+ω*V
    coeff = [DirectEvalLHSP10at0(n+2);DirectEvalLHSP10atPrime0(n+2)*Derivative(Jacobi(0,1,0..1),1)[1:n+2,1:n+2];V[1:n,1:n+2]] \ [0.0;0.0;pad(gF.coefficients,n)]
    return Fun(Jacobi(0,1, 0..1), coeff)
end

####
## Compute the solution with same parameters as discussed by Hale and in our paper
## This correponds to computing the solution in Figure 7(a).
ν = 3
μ = 2
ω = 20
u = solveSec54(2000,g,Kfun,ν,μ,ω)
plot(u,color=:black, legend=:false, xlabel = "x", ylabel = "u(x)", legendfontsize=12, tickfontsize=10, thickness_scaling = 1 , grid=:none)

####
## A basic loop to compute errors between the above high order solution and increasingly accurate approximations
function errorvals(max,step,u)
    errorvec = []
    for n=10:step:max
        diff = maximum(abs.([maximum(solveSec54(n,g,Kfun,ν,μ,ω)-u),minimum(solveSec54(n,g,Kfun,ν,μ,ω)-u)]))
        errorvec = push!(errorvec,diff)
    end
    return errorvec
end
####
## Plot successive errors to produce something like Figure 7(b). The exact shape of the plot depends on how dense we make the steps in n.
## Note that for low polynomial orders, since the order is too low to resolve the oscillations, the error varies a lot.
## Only once the polynomial order can resolve the oscillations do we see sensible convergence.
errorvec = errorvals(450,40,u)
plot(Array(10:40:450),errorvec,color=:black, legend=:bottomleft , yscale=:log10 , xlabel = "n", label="sparse method" , ylabel = "error", legendfontsize=12, tickfontsize=10, thickness_scaling = 1.2 , grid=:none)
