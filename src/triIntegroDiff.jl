#############################
## This file extends the linear functionality in src/triVolterraFullKernel.jl
## to nonlinear and integro-differential equations. This has many different variations.
## We hard-coded some standard ones here but it is straightforward to modify it for more general cases.

#############################
##   Provides a measure of distance between guess and true solution of non-linear 2nd kind Volterra where f(u) is the nonlinearity
####
function triNonLinearVolterraObjective(uguess,f,K,g::Fun,depth,kerneldepth)
    u = Fun(Jacobi(0,1,0..1),uguess)
    fofu = pad(Fun(x->f(u(x)),Jacobi(0,1,0..1)).coefficients,depth)
    V = Conversion(Jacobi(0,0,0..1),Jacobi(0,1,0..1))[1:depth,1:depth]*reflectPabtoPba(depth)*WLoweringP01P00(depth)*triVolterraFullKernelOpP01(K,depth,true,kerneldepth)
    return pad(uguess,depth)-V*fofu-pad(g.coefficients,depth)
end

#############################
##   Provides a measure of distance between guess and true solution of non-linear, integro-differential
##   2nd kind Volterra where f(u) is the nonlinearity, DD is the order of the derivative and
##   Init is a vector of length D containing the initial conditions.
##
##   Currently hard-coded for 1st and 2nd derivative problems only but easy to extend to general case when needed.
####
function triNonLinearIntegroDiffVolterraObjective(uguess,f,K,g,depth,DD,init,lincomb,kerneldepth,constkernel::Bool)
    u = Fun(Jacobi(0,1,0..1),uguess)
    gF = Fun(x->g(x),Jacobi(DD,DD+1, 0..1))
    fofu = pad(Fun(x->lincomb(x)+f(u(x)),Jacobi(0,1,0..1)).coefficients,depth)
    if constkernel # Check if kernel is a constant, in which case explicit monomial method is used
        sp = JacobiTriangle()
        xy = axes(sp,1)
        x,y = first.(xy),last.(xy)
        V = triVolterraOpP01([K(0.5,0.5),0.,0.,0.],depth,true)
    else
        V = triVolterraFullKernelOpP01(K,depth,true,kerneldepth)
    end
    V = reflectPabtoPba(depth)*WLoweringP01P00(depth)*V
    V = Conversion(Jacobi(0,0,0..1),Jacobi(DD,1+DD,0..1))[1:depth,1:depth]*V[1:depth,1:depth]
    if DD==1
        return pad([u(0)-init[1];pad((Derivative(DD)*u).coefficients,depth)-V*fofu-pad(gF.coefficients,depth)],depth)
    end
    if DD==2
        return pad([u(0)-init[1];(Derivative(DD-1)*u)(0)-init[2];pad((Derivative(DD)*u).coefficients,depth)-V*fofu-pad(gF.coefficients,depth)],depth)
    end
end
#############################
##   This is the pochhammer symbol function from FastTransforms.jl. It's reproduced here directly for testing purposes.
####
function pochhammer(x::Number,n::Integer)
    ret = one(x)
    if n≥0
        for i=0:n-1
            ret *= x+i
        end
    else
        ret /= pochhammer(x+n,-n)
    end
    ret
end

#############################
##   alternI: Diagonal matrix with entries (-1)^n
####
function alternI(N::Integer)
    altI=BandedMatrix{Float64}(undef, (N, N), (0,0))
    for n=1:N
        altI[n,n]=(-1)^(n+1)
    end
    return altI
end

#############################
##   EvalJPat1mat: Operator which evaluates a coefficient vector of
##              0..1 Jacobi polynomials in P_n^(α,β)basis at x=0 if given α.
##              Combine with alternI(N) and use β instead of α to evaluate at x=1.
##              Be mindful that in the ApproxFun version used here (though not in later versions), we have Jacobi(β,α) = P_n^(α,β)(x)
####
function EvalJPat1mat(α,N::Integer)
    poch=BandedMatrix{Float64}(undef, (N, N), (0,0))
    for n=1:N
        poch[n,n]=FastTransforms.pochhammer(1+α,n-1)/SpecialFunctions.gamma(n);
    end
    return poch
end
function EvalatLHS(β,N)
    return Ones(1,N)*alternI(N)*EvalJPat1mat(β,N)
end
function EvalatRHS(α,N)
    return Ones(1,N)*EvalJPat1mat(α,N)
end

#############################
## Some specific commonly occuring speccial cases of the above for higher efficiency computations, evaluate LHS of specific polynomial bases
## Prime in name indicates number of derivatives. Can also just use the general function but this is a bit faster if you know what you will need.
####
function DirectEvalLHSP10at0(N::Integer)
    vek=ones(1,N)
    for n=1:N
        vek[1,n]=(-1)^(n-1)
    end
    return vek
end
function DirectEvalLHSP10atPrime0(N::Integer)
    vek=ones(1,N)
    for n=1:N
        vek[1,n]=(-1)^(n-1)*n
    end
    return vek
end
function DirectEvalLHSP10atPPrime0(N::Integer)
    vek=ones(1,N)
    for n=1:N
        vek[1,n]=(-1)^(n-1)*n*(n+1)/2
    end
    return vek
end
function DirectEvalLHSP10atPPPrime0(N::Integer)
    vek=ones(1,N)
    for n=1:N
        vek[1,n]=(-1)^(n-1)*n*(n+1)*(n+2)/6
    end
    return vek
end
function DirectEvalLHSP10atPPPPrime0(N::Integer)
    vek=ones(1,N)
    for n=1:N
        vek[1,n]=(-1)^(n-1)*(n+3)*(n+2)*(n+1)*n/24
    end
    return vek
end
