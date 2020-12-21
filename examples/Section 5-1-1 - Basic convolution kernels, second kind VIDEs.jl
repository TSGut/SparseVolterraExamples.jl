#############################
##  Numerical experiments associated with section 5.1.1
#############################

using ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, SpecialFunctions, Plots
using SparseVolterraExamples

#############################
## Problem in Equation (11)
#############################
####
## This function sets up g in the correct basis as well as kernel K.
function solveSec511A(n)
    gf(x) = 1
    gF = Fun(x->gf(x),Jacobi(2,3, 0..1))
    Kfun(x,y) = x-y
        V = triVolterraFullKernelOpP01(Kfun,n,true)
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V
        V = Conversion(Jacobi(0,0,0..1),Jacobi(2,3,0..1))[1:n,1:n]*V[1:n,1:n]
        V = Derivative(Jacobi(0,1,0..1),2)[1:n,1:n]-V
    coeff = [DirectEvalLHSP10at0(n);DirectEvalLHSP10atPrime0(n)*Derivative(Jacobi(0,1,0..1),1)[1:n,1:n];V[1:n,1:n]] \ [1;0;pad(gF.coefficients,n)]
end

####
## Compute the solution
coeff = solveSec511A(50)
u = Fun(Jacobi(0,1,0..1), coeff)
####
## Plot the computed and analytic solution
plot(u,grid=false,ylabel="u(x)",xlabel="x",label="sparse method",legend=:topleft)
plot!(x->cosh(x),0,1,grid=false,ylabel="u(x)",xlabel="x",label="analytic cosh(x)",legend=:topleft)
####
## Plot the error of the computed compared to analytic solution.
plot(x->abs(cosh(x)-u(x)),0,1,grid=false,ylabel="error",xlabel="x")
####
## Looping over the the degree of polynomials gives a convergence plot like in Figure 1.
function convSec511A(max,stepsize)
    errorvec = []
    for n=1:stepsize:max
        diff = maximum(abs.([maximum(Fun(Jacobi(0,1,0..1),solveSec511A(n))-Fun(x->cosh(x),0..1)),minimum(Fun(Jacobi(0,1,0..1),solveSec511A(n))-Fun(x->cosh(x),0..1))]))
        errorvec = append!(errorvec,diff)
    end
    return errorvec
end
plot(Array(1:1:14),convSec511A(14,1),grid=false,ylabel="error",xlabel="n",yaxis=:log10,label=false)

#############################
## Problem in Equation (12)
#############################
####
## This function sets up g in the correct basis as well as kernel K.
function solveSec511B(n)
    gf(x) = -1+x
    gF = Fun(x->gf(x),Jacobi(4,5, 0..1))
    Kfun(x,y) = y-x
    V = triVolterraFullKernelOpP01(Kfun,n,true)
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V
        V = Conversion(Jacobi(0,0,0..1),Jacobi(4,5,0..1))[1:n,1:n]*V[1:n,1:n]
        V = Derivative(Jacobi(0,1,0..1),4)[1:n,1:n]-V
    coeff = [DirectEvalLHSP10at0(n);DirectEvalLHSP10atPrime0(n)*Derivative(Jacobi(0,1,0..1),1)[1:n,1:n];DirectEvalLHSP10atPPrime0(n)*Derivative(Jacobi(0,1,0..1),2)[1:n,1:n];DirectEvalLHSP10atPPPrime0(n)*Derivative(Jacobi(0,1,0..1),3)[1:n,1:n];V[1:n,1:n]] \ [-1;1;1;-1;pad(gF.coefficients,n)]
end
####
## Compute the solution
coeff = solveSec511B(15)
u = Fun(Jacobi(0,1,0..1), coeff)
####
## Plot the computed and analytic solution
plot(u,grid=false,ylabel="u(x)",xlabel="x",label="sparse method",legend=:topleft)
plot!(x->sin(x)-cos(x),0,1,grid=false,ylabel="u(x)",xlabel="x",label="analytic sin(x)-cos(x)",legend=:topleft)
####
## Plot the error of the computed compared to analytic solution.
plot(x->abs(sin(x)-cos(x)-u(x)),0,1,grid=false,ylabel="error",xlabel="x")
####
## Looping over the the degree of polynomials gives a convergence plot like in Figure 1.
function convSec511B(max,stepsize)
    errorvec = []
    for n=1:stepsize:max
        diff = maximum(abs.([maximum(Fun(Jacobi(0,1,0..1),solveSec511B(n))-Fun(x->sin(x)-cos(x),0..1)),minimum(Fun(Jacobi(0,1,0..1),solveSec511B(n))-Fun(x->sin(x)-cos(x),0..1))]))
        errorvec = append!(errorvec,diff)
    end
    return errorvec
end
plot(Array(1:4:40),convSec511B(40,4),grid=false,ylabel="error",xlabel="n",yaxis=:log10,label=false)

#############################
## Problem in Equation (13)
#############################
####
## This function sets up g in the correct basis as well as kernel K.
function solveSec511C(n)
    gf(x) = 1+x+(x^2)/2-(x^4)/(4*3*2)
    gF = Fun(x->gf(x),Jacobi(3,4, 0..1))
    Kfun(x,y) = (x-y)^2/2
    V = triVolterraFullKernelOpP01(Kfun,n,true)
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V
        V = Conversion(Jacobi(0,0,0..1),Jacobi(3,4,0..1))[1:n,1:n]*V[1:n,1:n]
        V = Derivative(Jacobi(0,1,0..1),3)[1:n,1:n]-V
    coeff = [DirectEvalLHSP10at0(n);DirectEvalLHSP10atPrime0(n)*Derivative(Jacobi(0,1,0..1),1)[1:n,1:n];DirectEvalLHSP10atPPrime0(n)*Derivative(Jacobi(0,1,0..1),2)[1:n,1:n];V[1:n,1:n]] \ [1;2;1;pad(gF.coefficients,n)]
end
####
## Compute the solution
coeff = solveSec511C(15)
u = Fun(Jacobi(0,1,0..1), coeff)
####
## Plot the computed and analytic solution
plot(u,grid=false,ylabel="u(x)",xlabel="x",label="sparse method",legend=:topleft)
plot!(x->x+exp(x),0,1,grid=false,ylabel="u(x)",xlabel="x",label="analytic x+exp(x)",legend=:topleft)
####
## Plot the error of the computed compared to analytic solution.
plot(x->abs(x+exp(x)-u(x)),0,1,grid=false,ylabel="error",xlabel="x")
####
## Looping over the the degree of polynomials gives a convergence plot like in Figure 1.
function convSec511C(max,stepsize)
    errorvec = []
    for n=1:stepsize:max
        diff = maximum(abs.([maximum(Fun(Jacobi(0,1,0..1),solveSec511C(n))-Fun(x->x+exp(x),0..1)),minimum(Fun(Jacobi(0,1,0..1),solveSec511C(n))-Fun(x->x+exp(x),0..1))]))
        errorvec = append!(errorvec,diff)
    end
    return errorvec
end
plot(Array(1:1:20),convSec511C(20,1),grid=false,ylabel="error",xlabel="n",yaxis=:log10,label=false)
