using ApproxFun, BandedMatrices, BlockBandedMatrices, LinearAlgebra, BlockArrays, MultivariateOrthogonalPolynomials, Test
using SparseVolterraExamples

#############################
## Tests for Volterra integral equations of first and second kind
## with Clenshaw Jacobi expansion of the kernel.
####
## Includes tests for 1 to x and 1 to 1-x limits of integration.
#############################

@testset "Volterra Integral Equation 1st kind, could be monomial, flip=true" begin
    g(x) = sin(x)-x*cos(x);
    K(x,y) = y;
    @test triVolterraEQ1FullKernelSolver(g,K,30,true) ≈ Fun(x->sin(x),Jacobi(0,1, 0..1))
    #########################################################################
    g(x) = 2+x-2*exp(x)+x*exp(x)
    K(x,y) = x-y;
    @test triVolterraEQ1FullKernelSolver(g,K,20,true) ≈ Fun(x->x*exp(x),Jacobi(0,1, 0..1))
end

@testset "Volterra Integral Equation 1st kind, could be monomial, flip=false" begin
    N = 20
    g(x) = 3-x-2*exp(1-x)+(1-x)*exp(1-x);
    K(x,y) = 1-y-x;
    @test triVolterraEQ1FullKernelSolver(g,K,N,false) ≈ Fun(x -> x*exp(x), Jacobi(0,1,0..1))
end

@testset "Volterra Integral Equation 1st kind, not monomial, flip=false" begin
    N = 20;
    g(x) = 1/12*(4-4*(-1+x)^3+3*(-1+x)^4-4*cos((1-x)^3))
    K(x,y) = 1+y+sin(y^3);
    @test triVolterraEQ1FullKernelSolver(g,K,N,false) ≈ Fun(x->x^2,Jacobi(0,1,0..1))
end

@testset "Volterra Integral Equation 2nd kind, could be monomial, flip=false" begin
    N = 15;
    S = Jacobi(0,1,0..1);
    g = Fun(x -> (x-1)*cos(2*pi*(1-x))/(2*pi)+((4*pi^2+1)/(4*pi^2))*sin(2*pi*(1-x)), S);
    K(x,y)=y;
    @test triVolterraEQ2FullKernelSolver(g,K,N,false) ≈ Fun(x -> sin(2*pi*(1-x)), S)
    #########################################################################
    N = 15;
    S = Jacobi(0,1,0..1);
    g = Fun(x -> (-17/30)*x^6+2*x^5-(3/2)*x^4-(11/2)*x^3+(15/2)*x^2-(21/10)*x+1/6, S);
    K(x,y) = 7*x*y+5*y^2;
    @test triVolterraEQ2FullKernelSolver(g,K,N,false) ≈ Fun(x -> 1-x^3, S)
end

@testset "Volterra Integral Equation 2nd kind, not monomial, flip=true" begin
    g = Fun(x->x^2+2*x-exp(x)+2, Jacobi(0,1,0..1));
    K(x,y) = (y-x)^2;
    @test triVolterraEQ2FullKernelSolver(g,K,20,true) ≈ Fun(x->exp(x),Jacobi(0,1,0..1));
end
