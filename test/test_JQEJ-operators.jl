using ApproxFun, BandedMatrices, BlockBandedMatrices, SparseArrays, LinearAlgebra, LazyArrays, BlockArrays, MultivariateOrthogonalPolynomials, Test
using SparseVolterraExamples

#############################
## Various tests for integration Qy and expansion Ey operator consistency.
#############################
@testset "Testing that E_x from P01 works correctly" begin
    g = Fun(x -> exp(x+3x^2), Jacobi(0,1,0..1));
    N = ncoefficients(g);
    G = Fun(JacobiTriangle(), triExgenP01(N) * g.coefficients)
    @test G(0.1,0.2) ≈ g(0.1);
    #########################################################################
    f = Fun(x -> exp(x+3x^2), Legendre(0..1));
    N = ncoefficients(f);
    C = (I : Legendre(0..1) → Jacobi(0,1,0..1))[1:N,1:N];
    F = Fun(JacobiTriangle(), triExgenP01(N) * (C * f.coefficients));
    @test F(0.1,0.2) ≈ f(0.1);
end

@testset "Testing that Q and E_y from P01 work correctly" begin
    # define a function which we want to integrate from 0 to 0.9
    # the expected result for this function is approx. 5.24914
    f = Fun(y -> exp(y+3*y^2), Jacobi(0,1,0..1));
    Ey=triEygenP01(length(f.coefficients));
    F = Fun(JacobiTriangle(), Ey*f.coefficients);
    @test F(0.1,0.2) ≈ f(0.2);
    @test F(0.1,0.3) ≈ f(0.3);
    #########################################################################
    # define a function which we want to integrate from 0 to 0.9
    # the expected result for this function is approx. 5.24914
    f = Fun(y -> exp(y+3*y^2), Jacobi(0,1,0..1));
    N = ncoefficients(f);
    # option 1: use direct ApproxFun + Q generator. This does not use E.
    @test triNoKernelDirect(f,0.9) ≈ sum(Fun(x -> f(x), 0..0.9));
    # option 2: generate Q and E separately and apply them to f.coefficients
    Q=triQgen(N);
    Ey=triEygenP01(length(f.coefficients));
    @test Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*f.coefficients)(0.1) ≈
                sum(Fun(x -> f(x), 0..0.9));
    # we can check that Ey does what it should like this:
    @test pad(Fun((x,y) -> f(y), JacobiTriangle()).coefficients,length(Ey*f.coefficients)) ≈ Ey*f.coefficients;
    #########################################################################
    # again for a different function from 0 to 0.5, expected result here is approx. 0.0917462
    g = Fun(y -> sin(8*y)+y^3-tan(y), Jacobi(0,1,0..1));
    N = ncoefficients(g);

    # option 1: use direct ApproxFun + Q generator. This does not use E.
    triNoKernelDirect(g,0.5);

    # option 2: generate Q and E separately and apply them to f.coefficients
    Q=triQgen(N);
    Ey=triEygenP01(N);
    @test Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*g.coefficients)(0.5) ≈
                sum(Fun(x -> g(x), 0..0.5));

    # we can check that Ey does what it should like this:
    @test pad(Fun((x,y) -> g(y), JacobiTriangle()).coefficients,length(Ey*g.coefficients)) ≈ Ey*g.coefficients;
end

@testset "Testing that Q and E_y from P00 work correctly" begin
    # define a function which we want to integrate from 0 to 0.9
    # the expected result for this function is approx. 5.24914
    f = Fun(y -> exp(y+3*y^2), Legendre(0..1));
    Ey=triEygenP00(length(f.coefficients));
    F = Fun(JacobiTriangle(), Ey*f.coefficients);
    @test F(0.1,0.2) ≈ f(0.2);
    @test F(0.1,0.3) ≈ f(0.3);
    #########################################################################
    # define a function which we want to integrate from 0 to 0.9
    # the expected result for this function is approx. 5.24914
    f = Fun(y -> exp(y+3*y^2), Legendre(0..1));
    N = ncoefficients(f);
    # option 1: use direct ApproxFun + Q generator. This does not use R.
    @test triNoKernelDirect(f,0.9) ≈ sum(Fun(x -> f(x), 0..0.9));
    # option 2: generate Q and E separately and apply them to f.coefficients
    Q=triQgen(N);
    Ey=triEygenP00(length(f.coefficients));
    @test Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*f.coefficients)(0.1) ≈
                sum(Fun(x -> f(x), 0..0.9));
    # we can check that Ey does what it should like this:
    @test pad(Fun((x,y) -> f(y), JacobiTriangle()).coefficients,length(Ey*f.coefficients)) ≈ Ey*f.coefficients;
    #########################################################################
    # again for a different function from 0 to 0.5, expected result here is approx. 0.0917462
    g = Fun(y -> sin(8*y)+y^3-tan(y), Legendre(0..1));
    N = ncoefficients(g);

    # option 1: use direct ApproxFun + Q generator. This does not use R.
    triNoKernelDirect(g,0.5);

    # option 2: generate Q and R separately and apply them to f.coefficients
    Q=triQgen(N);
    Ey=triEygenP00(N);
    @test Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*g.coefficients)(0.5) ≈
                sum(Fun(x -> g(x), 0..0.5));

    # we can check that Ey does what it should like this:
    @test pad(Fun((x,y) -> g(y), JacobiTriangle()).coefficients,length(Ey*g.coefficients)) ≈ Ey*g.coefficients;
end

#############################
## Various tests for the interactions of Qy,Ey with the Jacobi operators
#############################

@testset "Testing that Q*E*f integrates correctly" begin
    S = Jacobi(0,1,0..1);
    f = Fun(x -> 1, S);
    N = length(f.coefficients);
    Q = triQgen(N);
    Ey = triEygenP01(N);
    QEf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*f.coefficients);
    @test triNoKernelDirect(f,1.) ≈ QEf(0);
    #########################################################################
    S = Jacobi(0,1,0..1);
    f = Fun(x -> x, S);
    N = length(f.coefficients);
    Q = triQgen(N);
    Ey = triEygenP01(N);
    QEf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*f.coefficients);
    @test triNoKernelDirect(f,1.) ≈ QEf(0);
    #########################################################################
    S = Jacobi(0,1,0..1);
    f = Fun(x -> 1,);
    xf = Fun(x -> x, S);
    N = length(f.coefficients);
    y = Fun(identity, S);
    Jy = Multiplication(y, S);
    Q = triQgen(N);
    Ey = triEygenP01(N);
    QEf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*f.coefficients);
    QExf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*Jy[Block.(1:N), Block.(1:N)]*f.coefficients);
    @test triNoKernelDirect(f,1.) ≈ QEf(0);
    #########################################################################
    S = Jacobi(0,1,0..1);
    f = Fun(x -> 1+exp(3*x^2), S);
    xf = Fun(x -> x*(1+exp(3*x^2)),S);
    N = length(f.coefficients);
    y = Fun(identity, S);
    Jy = Multiplication(y, S);
    Q = triQgen(N);
    Ey = triEygenP01(N);
    QEf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*f.coefficients);
    QExf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)),(1/2)*Q*Ey*Jy[Block.(1:N), Block.(1:N)]*f.coefficients);
    @test triNoKernelDirect(f,1.) ≈ QEf(0);
end

@testset "Testing if expected commutation relations hold" begin
    # we test if the direct Q*EyP01 diagonal matrix is computed correctly
    N=50;
    @test triQEygenP01(N) ≈ triQgen(N)*triEygenP01(N);
    #########################################################################
    # we test expected commutation relations of R and Y on the triangle
    # once with a function and one without function
    N = 50;
    S = Jacobi(0,1,0..1);
    T = JacobiTriangle();
    f = Fun(x->2+x+exp(x^5+1), S);
    E = sparse(triEygenP01(N));
    y = Fun(identity, S);
    J = sparse(Multiplication(y, S)[Block.(1:N), Block.(1:N)]);
    x, y = Fun(identity, T);
    Y = sparse(Multiplication(y, T)[Block.(1:N), Block.(1:N)]);
    @test Y*E ≈ E*J;
    @test Y*E*pad(f.coefficients,N) ≈ E*(J*pad(f.coefficients,N));
    #########################################################################
    # we test expected commutation relations of Q and X
    N = 50;
    S = Jacobi(0,1,0..1);
    T = JacobiTriangle();
    Q = sparse(triQgen(N));
    x = Fun(identity, S);
    J = sparse(Multiplication(x, S)[Block.(1:N), Block.(1:N)]);
    x, y = Fun(identity, T);
    X = sparse(Multiplication(x, T)[Block.(1:N), Block.(1:N)]);
    @test Q*X ≈ J*Q
end

@testset "Selective testing of monomial kernel basis consistency" begin
    k = 81;
    N = 30;
    basisvec = triJQEJbasislexP01(k,N,false);
    basisgrid = triJQEJbasisgridgenP01(k,k,N);
    # first couple
    @test basisvec[1] ≈ basisgrid[1][1];
    @test basisvec[2] ≈ basisgrid[1][2];
    @test basisvec[3] ≈ basisgrid[2][1];
    @test basisvec[4] ≈ basisgrid[2][2];
    # square number sequence elements
    @test basisvec[9] ≈ basisgrid[3][3];
    @test basisvec[16] ≈ basisgrid[4][4];
    @test basisvec[25] ≈ basisgrid[5][5];
    @test basisvec[36] ≈ basisgrid[6][6];
    @test basisvec[81] ≈ basisgrid[9][9];
    # some pseudo random picks
    @test basisvec[12] ≈ basisgrid[3][4];
    @test basisvec[27] ≈ basisgrid[2][6];
    @test basisvec[80] ≈ basisgrid[9][1];
end
