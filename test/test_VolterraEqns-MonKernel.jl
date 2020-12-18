using ApproxFun, BandedMatrices, BlockBandedMatrices, LinearAlgebra, BlockArrays, MultivariateOrthogonalPolynomials, Test
using SparseVolterraExamples

#############################
## Tests for Volterra integral equations of first and second kind
## with non-Clenshaw monomial expansion of the kernel.
####
## Includes tests for 1 to x and 1 to 1-x limits of integration.
#############################

@testset "Volterra Integral Equation 1st kind, mon. kernels, flip=false" begin
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> 1-sin(x)^2, S)
    # K given in monomial basis
    K = [2.3,0.33,2.1,1.2,7.6]  # = K(x,y)=2.3+0.33*y+2.1*x+1.2*x*y+7.6*y^2
        V=triVolterraOpP01(K,N,false)
        x = Fun(identity,S)
        J = Multiplication(x, S)[1:N,1:N]
    VintFun1=Fun(Jacobi(0,1,0..1),(V-J*V)*pad(g.coefficients,N))
    VintFun2=triVolterraFunP01(g,K,N,false)
    Vtest=triVolterraEQ1Solver(VintFun1,K,N,false)
    @test VintFun1 ≈ VintFun2
    @test Vtest ≈ g
    #########################################################################
    # sin(1-x) = integral_0^(1-x) u(t) dt
    # analytic solution is u(x)=cos(x)
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> sin(1-x), S)
    K = [1.,0.,0.,0.]  # = K(x,y)=1
    @test triVolterraEQ1Solver(g,K,N,false) ≈ Fun(x -> cos(x), S)
    #########################################################################
    # 2-x-exp(1-x) = integral_0^(1-x) (y-1-x) u(y) dy
    # analytic solution is u(x)=exp(x)
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> 2-x-exp(1-x), S)
    K = [-1.0,1.0,1.0,0.]  # = K(x,y)=(y-1+x)
    @test triVolterraEQ1Solver(g,K,N,false) ≈ Fun(x -> exp(x), S)
    #########################################################################
    # 3-x-2*exp(1-x)+(1-x)*exp(1-x) = integral_0^(1-x) (1-x-y) u(y) dy
    # analytic solution is u(x)=x*exp(x)
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> 3-x-2*exp(1-x)+(1-x)*exp(1-x), S)
    K = [1.,-1.,-1.0,0.]  # = K(x,y)=(1-x-y)
    @test triVolterraEQ1Solver(g,K,N,false) ≈ Fun(x -> x*exp(x), S)
end

@testset "Volterra Integral Equation 1st kind, mon. kernels, flip=true" begin
    g=Fun(x -> sin(x)-(x)*cos(x), Jacobi(0,1,0..1))
    K = [0.,1.,0.,0.]
    N = 20
    @test triVolterraEQ1Solver(g,K,N,true)  ≈ Fun(x->sin(x),Jacobi(0,1,0..1))
    #########################################################################
    g=Fun(x -> 2+x-2*exp(x)+x*exp(x), Jacobi(0,1,0..1))
    K = [0.,-1.,1.,0.]
    N = 15
    @test triVolterraEQ1Solver(g,K,N,true) ≈ Fun(x->x*exp(x),Jacobi(0,1,0..1))
end

@testset "Volterra Integral Equation 2nd kind, mon. kernels, flip=false" begin
    # u(x) = (x-1)*cos(2*pi*(1-x))/(2*pi)+((4*pi^2+1)/(4*pi^2))*sin(2*pi*(1-x))+\int_0^(1-x) t*u(t)dt
    # solution is u(x) = sin(2*pi*(1-x))
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> (x-1)*cos(2*pi*(1-x))/(2*pi)+((4*pi^2+1)/(4*pi^2))*sin(2*pi*(1-x)), S)
    K = zeros(5)
    K[2]=1.
    @test triVolterraEQ2Solver(g,K,N,false) ≈ Fun(x -> sin(2*pi*(1-x)), S)
    #########################################################################
    # u(x) = (-17/30)*x^6+2*x^5-(3/2)*x^4-(11/2)*x^3+(15/2)*x^2-(21/10)*x+1/6 + \int_0^(1-x) (7xy+5y^2) u(y)dy
    # solution is u(x) = 1-x^3
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> (-17/30)*x^6+2*x^5-(3/2)*x^4-(11/2)*x^3+(15/2)*x^2-(21/10)*x+1/6, S)
    K = [0.,0.,0.,7.,5.,0.,0.,0.,0.]
    @test triVolterraEQ2Solver(g,K,N,false) ≈ Fun(x -> 1-x^3, S)
    #########################################################################
    # u(x) = -x+x^2/2+1/2 + \int_0^(1-x) u(t)dt
    # solution is u(x) = 1-x
    N=15
    S = Jacobi(0,1,0..1)
    g = Fun(x -> -x+x^2/2+1/2, S)
    K = [1.,0.,0.]
    @test triVolterraEQ2Solver(g,K,N,false) ≈ Fun(x -> 1-x, S)
end

@testset "Volterra Integral Equation 2nd kind, mon. kernels, flip=true" begin
    g=Fun(x -> 1, Jacobi(0,1,0..1))
    K = [0.,1.,-1.,0.]
    N = 10
    @test triVolterraEQ2Solver(g,K,N,true)  ≈ Fun(x->cos(x),Jacobi(0,1,0..1))
    #########################################################################
    g=Fun(x -> 1-x-0.5*x^2, Jacobi(0,1,0..1))
    K = [0.,-1.,1.,0.]
    N = 10
    @test triVolterraEQ2Solver(g,K,N,true)  ≈ Fun(x->1-sinh(x),Jacobi(0,1,0..1))
    #########################################################################
    g=Fun(x -> x+x^4+0.5*x^2+0.2*x^5, Jacobi(0,1,0..1))
    K = [-1.,0.,0.,0.]
    N = 25
    @test triVolterraEQ2Solver(g,K,N,true)  ≈ Fun(x->x+x^4,Jacobi(0,1,0..1))
end
