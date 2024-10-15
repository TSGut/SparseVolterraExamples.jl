#############################
## The functions in this file may currently not work for some easy edge cases like K=constant.
## Use the monomial implementation in src/triVolterraMonomialKernel.jl for those cases instead
## and reserve this for non-trivial kernels.
#############################
##  triVolterraEQ1FullKernelSolver : Solves Volterra integral equation of the first kind, i.e.
##                         g(x)= \int_0^(x) K(x,t) u(t) dt
##                                          OR
##                         g(x) =  \int_0^(1-x) K(x,t) u(t) dt (if flip = false)
##                         where g and the kernel K are given and u(t) is the unknown
##                         to be solved for.
####
function triVolterraEQ1FullKernelSolver(g::Function,Kfun::Function,depth::Integer,flip::Bool,kerneldepth::Int)
    V = WLoweringP01P00(depth)*triVolterraFullKernelOpP01(Kfun,depth,flip,kerneldepth)
    x = Fun(identity,Jacobi(0,1,0..1))
    if flip
        gF = Fun(x->g(1-x),Jacobi(0,0, 0..1))
        u = Fun(Jacobi(0,1, 0..1), V\pad(gF.coefficients,depth))
    else
        gF = Fun(x->g(x),Jacobi(0,0, 0..1))
        u = Fun(Jacobi(0,1, 0..1), V\pad(gF.coefficients,depth))
    end
    return u
end
#############################
##  triVolterraEQ2FullKernelSolver : Solves Volterra integral equation of the second kind, i.e.
##                         u(x) = g(x) + \int_0^(x) K(x,t) u(t) dt
##                                          OR
##                         u(x) = g(x) + \int_0^(1-x) K(x,t) u(t) dt (if flip = false)
##                         where g and the kernel K are given and u(x) is the unknown
##                         to be solved for.
##
##                         This solver uses a weighted lowering operator approach (I-SRLV).
####
function triVolterraEQ2FullKernelSolver(g::Fun,Kfun::Function,depth::Integer,flip::Bool,kerneldepth::Int)
    V = triVolterraFullKernelOpP01(Kfun,depth,flip,kerneldepth)
    if flip
        V = I-Conversion(Jacobi(0,0,0..1),Jacobi(0,1,0..1))[1:depth,1:depth]*reflectPabtoPba(depth)*WLoweringP01P00(depth)*V
        u = Fun(Jacobi(0,1, 0..1), V\pad(g.coefficients,depth))
    else
        x = Fun(identity,Jacobi(0,1,0..1))
        V = I-V+Multiplication(x, Jacobi(0,1,0..1))[1:depth,1:depth]*V
        u = Fun(Jacobi(0,1, 0..1), V\pad(g.coefficients,depth))
    end
    return u
end
#############################
##   triVolterraFun  :  Returns a Fun which evaluated at x gives the Volterra integral
##                      of the function f given as a Fun.
##
##                      use "flip = false" for integral from 0 to 1-x
####
function triVolterraFullKernelFunP01(f::Fun,Kfun::Function,depth::Integer,flip::Bool,kerneldepth::Int)
        op = triVolterraFullKernelOpP01(Kfun,depth,flip,kerneldepth)
        S = Jacobi(0,1,0..1)
        x = Fun(identity,S)
        Jx = Multiplication(x, S)[1:depth,1:depth]
        cfs = (I-Jx)*op*pad(f.coefficients,depth)
        if flip
                return Fun(Jacobi(1,0,0..1),reflectPabtoPba(length(cfs))*cfs)
        else
                return Fun(Jacobi(0,1,0..1),cfs)
        end
end
#############################
##   triVolterraFullKernelOpP01:  This version allows explicit kernel depth input which sometimes improves performance
####
function triVolterraFullKernelOpP01(K::Function,depth::Integer,flip::Bool,kerneldepth::Integer)
        S = Jacobi(0,1,0..1)
        x = Fun(identity,S)
        Jxy = Multiplication(x, S)[1:depth,1:depth]
        sp = JacobiTriangle()
        xy = axes(sp,1)
        x,y = first.(xy),last.(xy)
        if flip
            Kfun = sp[:,Block.(1:kerneldepth)] \ (@. K(1 -x, y))
        else
            Kfun = sp[:,Block.(1:kerneldepth)] \ (@. K(x, y))
        end
        N = blocksize(Kfun)[1]
        cfs = BlockedArray(pad(Kfun, sum(1:N)),1:N)
        C = ClenshawRecurrenceData(sp,N+2)
        D = triQEygenP01(depth)
        B2 = Fill(D,N) .* view(cfs,Block(N))
        B1 = Fill(D,N-1) .* view(cfs,Block(N-1)) .+ Fill(Jxy,N-1) .* (C.B̃ˣ[N-1]*B2) .+ (C.B̃ʸ[N-1]*B2) .* Fill(Jxy,N-1) .- C.B[N-1]*B2
        for K = N-2:-1:1
            (B1, B2) =  Fill(D,K) .* view(cfs,Block(K)) .+ Fill(Jxy,K) .* (C.B̃ˣ[K]*B1) .+ (C.B̃ʸ[K]*B1) .* Fill(Jxy,K) .- C.B[K]*B1 .- C.C[K]*B2 ,  B1
        end
        return first(B1)
end
#############################################################
## The stuff below is currently not in use outside of testing but is supplied for completeness
#############################
#############################
##   triBandedSgen: generates the bounded and banded S^-1_(1,0)^(1,1)*reflectPabtoPba*S^_(1,0)^(1,1)*(1-J) directly
##                  as a tridiagonal BandedMatrix.
####
function triBandedSgen(N::Integer)
    QEy=BandedMatrix{Float64}(undef, (N, N), (1, 1))
    QEy[1,1] = 1/3
    QEy[2,1] = 1/3
    for c=2:N-1
        QEy[c-1,c] = (-1)^c * (c-1)/(4*c-2)                         #superdiagonal
        QEy[c,c] = (-1)^(c+1) * c/(4*c^2-1)                         #diagonal
        QEy[c+1,c] = (-1)^(c+1) * ((c-1)/(4*c-2)+c/(4*c^2-1))    #subdiagonal
    end
    QEy[N-1,N] = (-1)^N * (N-1)/(4*N-2)
    QEy[N,N] = (-1)^(N+1) * N/(4*N^2-1)
    return QEy
end
