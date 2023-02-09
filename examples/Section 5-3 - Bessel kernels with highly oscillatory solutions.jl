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
K(x,y) = besselj(μ,ω.*(x-y))
g(x) = besselj(μ+ν,ω*x)+1/(2*x^2)*((ν-1)*(ν-2)*besselj(ν-1,ω*x)+(ν+1)*(ν+2)*besselj(ν+1,ω*x))

#####
## Solver function with step-by-step explanation for given polynomial degree 'n'
## Note the factor 10^(-3) placed as discussed in the paper.
function solveSec54(n,g,Kfun,ν,μ,ω,M)
    gF = Fun(x->g(x),Jacobi(2,3, 0..1))
    V = triVolterraFullKernelOpP01(Kfun,n+2,true,M)
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
Kf(x,y) = besselj(2.,20. .*(x-y))
u = solveSec54(2000,g,Kf,ν,μ,ω,7)
plot(u,color=:black, legend=:false, xlabel = "x", ylabel = "u(x)", legendfontsize=12, tickfontsize=10, thickness_scaling = 1 , grid=:none)
