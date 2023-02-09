#############################
##  Numerical experiments associated with section 5.1.2
#############################
## While the paper does not explicitly treat third kind Volterra equations, an extension to these
## cases is relatively straightforward. Note that while in general there won't be exponential convergence
## unless the involved powers are favorably smooth, our method still achieves higher accuracy than many competitor methods.
##

using ApproxFun, MultivariateOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, SpecialFunctions, Plots
using SparseVolterraExamples

#############################
## Problem in Equation (17)
#############################
####
## The following block computes the coefficient vector of the approximation for a given polynomial order of approximation n.
## The accuracy obtained depends on both the polynomial order for the multiplication as well as the solution.
## Choosing the same order is fine in most cases but better accuracy can be obtained by adjusting for a given problem.
function solveSec512A(n,multin,M)
    gf(x) = x^(2/3)*(10/3*x^(7/3)-3/16*x^(14/3));
    Kfun(x,y) = y;
    V = triVolterraFullKernelOpP01(Kfun,n,true,M);
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V;
        V = Conversion(Jacobi(0,0,0..1),Jacobi(1,2,0..1))[1:n,1:n]*V;
        V = Multiplication(Fun(x->x^(2/3),Jacobi(1,2,0..1),multin),Jacobi(1,2,0..1))[1:n,1:n]*Derivative(Jacobi(0,1,0..1),1)[1:n,1:n]-V
    coeff = [DirectEvalLHSP10at0(n);V[1:n-1,1:n]] \ [0;pad(Fun(x->gf(x),Jacobi(1,2, 0..1)).coefficients,n-1)]
    return coeff
end
####
## Plot the analytic solution against the computed solution
coeff = solveSec512A(50,40,12)
plot(Fun(Jacobi(0,1,0..1),coeff),grid=false,xlabel="x",ylabel="u(x)",label="sparse method")
plot!(x->(x)^(10/3),0,1,grid=false,xlabel="x",ylabel="u(x)",label="analytic solution")
####
## Check the numerical error.
plot(x->(Fun(Jacobi(0,1,0..1),coeff)(x)-(x)^(10/3)),0,1,grid=false,xlabel="x",ylabel="error",label=false)

#############################
## Problem in Equation (18)
#############################
####
## The following block computes the coefficient vector of the approximation for a given polynomial order of approximation n.
## The accuracy obtained depends on both the polynomial order for the multiplication as well as the solution.
## Choosing the same order is fine in most cases but better accuracy can be obtained by adjusting for a given problem.
function solveSec512B(n,multin,M)
    gf(x) = 9/2*x^4-1/20*x^(11/2)-1/6*x^6;
    Kfun(x,y) = sqrt(y);
    V = triVolterraFullKernelOpP01(Kfun,n,true,M);
        V = reflectPabtoPba(n)*WLoweringP01P00(n)*V;
        V = Conversion(Jacobi(0,0,0..1),Jacobi(1,2,0..1))[1:n,1:n]*V;
        V = Multiplication(Fun(x->sqrt(x),Jacobi(1,2,0..1),multin),Jacobi(1,2,0..1))[1:n,1:n]*Derivative(Jacobi(0,1,0..1),1)[1:n,1:n]-(1/20)*Conversion(Jacobi(0,1,0..1),Jacobi(1,2,0..1))[1:n,1:n]*Multiplication(Fun(x->x,Jacobi(0,1,0..1),n),Jacobi(0,1,0..1))[1:n,1:n]-V
    coeff = [DirectEvalLHSP10at0(n);V[1:n-1,1:n]] \ [0;pad(Fun(x->gf(x),Jacobi(1,2, 0..1)).coefficients,n-1)]
    return coeff
end
####
## Plot the analytic solution against the computed solution
coeff = solveSec512B(80,30,13)
plot(Fun(Jacobi(0,1,0..1),coeff),grid=false,xlabel="x",ylabel="u(x)",label="sparse method")
plot!(x->(x)^(9/2),0,1,grid=false,xlabel="x",ylabel="u(x)",label="analytic solution")
####
## Check the numerical error.
plot(x->(Fun(Jacobi(0,1,0..1),coeff)(x)-(x)^(9/2)),0,1,grid=false,xlabel="x",ylabel="error",label=false)
