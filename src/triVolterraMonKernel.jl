#############################
## The functions in this file compute linear Volterra integral equations
## on the condition that the kernel can be supplied as a vector of polynomial coefficients.
## Its main use outside of constant kernel edge cases is testing.
## In most cases you will want to instead use the functions in src/triVolterraFullKernel.jl
#############################
##  triVolterraEQ1Solver : Solves Volterra integral equation of the first kind, i.e.
##                         g(x)= \int_0^(x) K(x,t) u(t) dt
##                                          OR
##                         g(x) =  \int_0^(1-x) K(x,t) u(t) dt (if flip = false)
##                         where g and the kernel K are given and u(t) is the unknown
##                         to be solved for. K must be monomial coefficient vector.
####
function triVolterraEQ1Solver(g,Kmon::Array{Float64,1},N::Integer,flip::Bool)
    V = triVolterraOpP01(Kmon,N,flip)
    x = Fun(identity,Jacobi(0,1,0..1))
    V = V-Multiplication(x, Jacobi(0,1,0..1))[1:N,1:N]*V
    if flip
        g = Fun(x->g(1-x),Jacobi(0,1, 0..1))
        u = Fun(Jacobi(0,1, 0..1), V\pad(g.coefficients,N))
    else
        u = Fun(Jacobi(0,1, 0..1), V\pad(g.coefficients,N))
    end
    return u
end
#############################
##  triVolterraEQ2Solver : Solves Volterra integral equation of the second kind, i.e.
##                         u(x) = g(x) + \int_0^(x) K(x,t) u(t) dt
##                                          OR
##                         u(x) = g(x) + \int_0^(1-x) K(x,t) u(t) dt (if flip = false)
##                         where g and the kernel K are given and u(x) is the unknown
##                         to be solved for. K must be monomial coefficient vector.
####
function triVolterraEQ2Solver(g,Kmon::Array{Float64,1},N::Integer,flip::Bool)
    V = triVolterraOpP01(Kmon,N,flip)
    x = Fun(identity,Jacobi(0,1,0..1))
    V = V-Multiplication(x, Jacobi(0,1,0..1))[1:N,1:N]*V
    if flip
        C01 = Conversion(Jacobi(0,1,0..1),Jacobi(1,1,0..1))[1:N,1:N]
        u = Fun(Jacobi(0,1, 0..1), (C01-reflectPabtoPba(N)*C01*V)\(C01*pad(g.coefficients,N)))
    else
        u = Fun(Jacobi(0,1, 0..1), (I-V)\pad(g.coefficients,N))
    end
    return u
end
#############################
##   triVolterraFun  :  Returns a Fun which evaluated at x gives the Volterra integral
##                      of the function f given as a Fun with the kernel K(x,y) given in monomial expansion.
##
##                      use "flip = false" for integral from 0 to 1-x
####
function triVolterraFunP01(f,Kmon::Array{Float64,1},N::Integer,flip::Bool)
    if length(Kmon)<2
        Kmon=pad(Kmon,2)
    end
    if flip
        x = Fun(identity,Jacobi(0,1,0..1))
        J = Multiplication(x, Jacobi(0,1,0..1))[1:N,1:N]
        VintFun = Fun(Jacobi(0,1, reverseorientation(0..1)), (I-J)*(triVolterraOpP01(Kmon,N,flip)*pad(f.coefficients,N)))
    else
        VintFun = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)), 1/2*triVolterraOpP01(Kmon,N,flip)*pad(f.coefficients,N))
    end
    return VintFun
end
#############################
##   triVolterraOpP01:  Generates Volterra integral operator sum(J^mQEJ^n)
##                      when given monomial expansion of kernel K(x,y)
##                      on JacobiTriangle()
####
function triVolterraOpP01(Kmon::Array{Float64,1},N::Integer,flip::Bool)
    # generating algorithm does not work for k=1 so we pad to k=2
    if length(Kmon)<2
        Kmon=pad(Kmon,2)
    end
    d_K = length(Kmon)
    k = Int(ceil(sqrt(d_K)))
    basis = triJQEJbasislexP01(k,N,flip)
    operator = Kmon[1]*basis[1]
    for j=2:d_K
        operator=operator+Kmon[j]*basis[j]
    end
    return operator
end
#############################
##   triJQEJbasislexP01:  Generates full basis of J^lQEJ^k in a lexicographically sorted vector.
##                        basis[1] will be 00-element, basis[2] will be 01 element and so on.
##                        k>1 is the highest order of x and y in the kernel.
####
function triJQEJbasislexP01(k::Integer,N::Integer,flip::Bool)
    # generating algorithm does not work for k=1 so we pad to k=2
    if k<2
        k=k+1
    end
    S = Jacobi(0,1,0..1)
    x = Fun(identity,S)
    Jy = Multiplication(x, S)[1:N,1:N]
    if flip
        Jx = Multiplication(1-x, S)[1:N,1:N]
    else
        Jx = Jy
    end
    QEy=triQEygenP01(N)
    basis=[]
    for j=1:(k+1)^2
        push!(basis,zeros(N,N))
    end
    # some initialization for convenience
    basis[1]=QEy
    basis[2]=QEy*Jy
    basis[3]=Jx*QEy
    for i=2:k
        # compute y-only
        basis[i^2+1]=basis[(i-1)^2+1]*Jy
        # compute x-only
        basis[(i+1)^2-1]=Jx*basis[i^2-1]
        # compute y-dominant
        for j=1:Int(i-1)
            basis[i^2+j+1]=Jx*basis[i^2+j]
        end
        # compute x-y-balanced
        basis[i^2]=Jx*basis[i^2-i]
        # compute x-dominant
        for j=1:Int(i-1)
            basis[(i+1)^2-j-1]=basis[(i+1)^2-j]*Jy
        end
    end
    basis[(k+1)^2]=Jx*basis[k^2-k]
    return basis
end
#############################################################
## The stuff below is currently not in use outside of testing
#############################
##   triJQEJbasisgridP01: Generates full basis of J^lQEJ^k in a grid.
##                        k and l are the highest orders of x and y in the kernel.
##                        basis[m][n] will give the m+1,n+1 element.
##                        Primarily for testing purposes, not recommended outside of testing.
####
function triJQEJbasisgridgenP01(k::Integer,l::Integer,N::Integer)
    S = Jacobi(0,1,0..1)
    x = Fun(identity,S)
    J = Multiplication(x, S)[1:N,1:N]
    QEy=triQEygenP01(N)
    # careful: basis[1][1] actually relects the 0,0 element and so on.
    basis=[]
    ybasisV = [QEy]
    for m=1:k
        for n=2:l
            push!(ybasisV,ybasisV[n-1]*J)
        end
        push!(basis,ybasisV)
        ybasisV=[J*basis[m][1]]
    end
    return basis
end
#############################
##   triNoKernelDirect: Direct integration from 0 to 1-x with K(x,y)=1 using only Q and ApproxFun, skipping the use of
##                      E and commutation relations.
##                      Primarily for testing purposes, not recommended outside of testing.
####
function triNoKernelDirect(f,x::Float64)
    # Computes integral of f(y) from 0 to 1-x via Jacobi polynomials on triangle.
    T = JacobiTriangle()
    f̃ = Fun((x,y) -> f(y), T)
        N = length(f̃.coefficients)
        v = PseudoBlockArray(pad(f̃.coefficients,sum(1:N)), 1:N)
    Q = triQgen(N)
    Qf = Fun(JacobiWeight(0,1, Jacobi(0,1, 0..1)), (1/2*Q*v))
    return Qf(1-x)
end
