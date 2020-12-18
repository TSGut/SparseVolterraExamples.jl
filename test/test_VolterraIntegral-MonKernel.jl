using ApproxFun, BandedMatrices, BlockBandedMatrices, LinearAlgebra, BlockArrays, MultivariateOrthogonalPolynomials, Test
using SparseVolterraExamples

#############################
## Tests for Volterra integral operator
## with non-Clenshaw monomial expansion of the kernel.
####
## Includes tests for 1 to x and 1 to 1-x limits of integration.
#############################

@testset "Volterra Integral operator testing, no kernel" begin
    N=15;
    S = Jacobi(0,1,0..1);
    g = Fun(x -> 1-sin(x)^2+x*exp(x^3), S);
    V = triVolterraFunP01(g,[1.,0.],N,false)(0);
    @test V ≈ triNoKernelDirect(g,1.);
end

@testset "Volterra Integral operator testing, with kernel" begin
    N=15;
    S = Jacobi(0,1,0..1);
    g = Fun(x -> 1-sin(x)^2, S);
    # K given in monomial basis
    K = [1.,0.,2.,1.]; # = K(x,y)=1+2*x+x*y
    val=0.4;
    Kfval = Fun(y -> 1+2*val+val*y, S);
    @test triVolterraFunP01(g,K,N,false)(val) ≈ triNoKernelDirect(Kfval*g,1.0-val);
    #########################################################################
    N=15;
    S = Jacobi(0,1,0..1);
    g = Fun(x -> 1-sin(x)^2, S);
    # K given in monomial basis
    K = [2.3,0.33,2.1,1.2,7.6]; # = K(x,y)=2.3+0.33*y+2.1*x+1.2*x*y+7.6*y^2
    val=0.6;
    Kfval = Fun(y -> 2.3+0.33*y+2.1*val+1.2*val*y+7.6*y^2, S);
    @test triVolterraFunP01(g,K,N,false)(val) ≈ triNoKernelDirect(Kfval*g,1.0-val);
end
