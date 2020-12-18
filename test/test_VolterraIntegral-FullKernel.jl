using ApproxFun, BandedMatrices, BlockBandedMatrices, LinearAlgebra, BlockArrays, MultivariateOrthogonalPolynomials, Test
using SparseVolterraExamples

#############################
## Tests for Volterra integral operator
## with Clenshaw Jacobi expansion of the kernel.
####
## Includes tests for 1 to x and 1 to 1-x limits of integration.
#############################

@testset "Set #1.1 - Volterra, Clenshaw sin kern, flip = true" begin
    K(x,y)=sin(x+y);
    f=Fun(y->y, Jacobi(0,1,0..1));
    @test triVolterraFullKernelFunP01(f,K,25,true) ≈ Fun(x->-x*cos(2*x)-sin(x)+sin(2*x), Jacobi(1,0,0..1));
end

@testset "Set #1.2 - Volterra, Clenshaw sin kern, flip = false, val test" begin
    K(x,y)=sin(x+y);
    val = rand(0..1);
    f=Fun(y->y, Jacobi(0,1,0..1));
    @test triVolterraFullKernelFunP01(f,K,30,false)(val) ≈ Fun(x->((-1)*cos(1)+x*cos(1)+sin(1)-sin(x)), Jacobi(1,0,0..1))(val);
end

@testset "Set #2 - Volterra, Clenshaw sinh kern, flip = true" begin
    K(x,y)=(y-x)*sinh(x);
    f=Fun(y->y, Jacobi(0,1,0..1));
    @test triVolterraFullKernelFunP01(f,K,25,true) ≈ Fun(x->-1/6*x^3*sinh(x), Jacobi(1,0,0..1));
end

@testset "Set #3 - Volterra, Clenshaw sinh kern, flip = false" begin
    K(x,y)=(y-x)*sinh(x^2);
    f=Fun(y->y, Jacobi(0,1,0..1));
    @test triVolterraFullKernelFunP01(f,K,25,false) ≈ Fun(x->(sinh(x^2))/3-3/2*x*sinh(x^2)+2*x^2*sinh(x^2)-5/6*x^3*sinh(x^2), Jacobi(0,1,0..1));
end
