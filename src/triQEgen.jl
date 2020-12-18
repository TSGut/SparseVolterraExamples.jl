#############################
##   triQgen: Generate Q operator which represents integration with regards to y
##               on JacobiTriangle().
####
function triQgen(N::Integer)
    # Generates Q operator as a block banded matrix.
    # Q*f(x,y) can be used to integrate f from 0 to 1-u via Jacobi polynomials on triangle.
    Q = BlockBandedMatrix(Zeros(N, sum(1:N)), Ones{Int}(N), 1:N, (0,0))
        for n=1:N
            view(Q, Block(n,n))[1,1] = 1.
        end
    return Q
end
#############################
##   triEygenP01: Extension operator for y from Jacobi P^(0,1) basis to canonical Jacobi triangle
####
function triEygenP01(N::Integer)
    Ey = BlockBandedMatrix(Zeros(sum(1:N),N), 1:N,Ones{Int}(N), (0,0))
        for n=1:N
            for j=1:n
                view(Ey, Block(n,n))[j,1] = (-1)^n*(-1)^j*(2*j-1)/n
            end
        end
    return Ey
end
#############################
##   triQEygenP01: This directly generates Q*EyP01 as a diagonal matrix
##                 without first computing the operators Q and RyP01
####
function triQEygenP01(N::Integer)
    QEy=BandedMatrix{Float64}(undef, (N, N), (0,0))
    for n=1:N
        QEy[n,n]=(-1)^(n+1)/n
    end
    return QEy
end
#############################
##   WLoweringP01P00: Weighted lowering operator which also multiplies by (1-x) from
##                      Jacobi(0,1) to Jacobi(0,0)
function WLoweringP01P00(N::Integer)
    WLow=I+BandedMatrix{Float64}(-0.5*Ones(N,N), (1,0))
    return WLow
end
#############################
##   reflectPabtoPba: Reflects coefficient vector from P^(a,b)(1-x) basis into P^(b,a)(x) basis
####
function reflectPabtoPba(N::Integer)
    Refl=BandedMatrix{Float64}(undef, (N, N), (0,0))
    for n=1:N
        Refl[n,n]=(-1)^(n+1)
    end
    return Refl
end
#############################################################
## The functions below are not in use outside of testing purposes
#############################
##   WLoweringP01Pm11: Weighted lowering operator which also multiplies by (1-x) from
##                      Jacobi(0,1) to Jacobi(1,-1)
function WLoweringP01Pm11(N::Integer)
    WLow=BandedMatrix{Float64}(undef, (N, N), (1,0))
    for c=1:N-1
        WLow[c,c]=(1-c)/(2*c)
        WLow[c+1,c]=-0.5
    end
    WLow[N,N]=(1-N)/(2*N)
    return WLow
end
##   ReflectorP01: Reflects coefficient vector from P^(1,0)(x) basis into P^(1,0)(1-t) basis
##      avoid use because it is upper triangular
####
function ReflectorP01(N::Integer)
    Reflector=BandedMatrix{Float64}(undef, (N, N), (0,N))
    for r=1:N
        for col=r+1:N
            Reflector[r,col]=(-1)^(r+1)*2*r/col
        end
            Reflector[r,r]=(-1)^(r+1)
    end
    return Reflector
end
#############################
##   triEygenP00: Extension operator for y from Legendre P^(0,0) basis to canonical Jacobi triangle
####
function triEygenP00(N::Integer)
    Ey = BlockBandedMatrix(Zeros(sum(1:N),N), 1:N,Ones{Int}(N), (0,1))
        for n=1:N
            for j=1:n-1
                view(Ey, Block(n-1,n))[j,1] = (-1)^n*(-1)^j*(2*j-1)/(2*n-1)
                view(Ey, Block(n,n))[j,1] = (-1)^n*(-1)^j*(2*j-1)/(2*n-1)
            end
            view(Ey, Block(n,n))[n,1] = 1.0
        end
    return Ey
end
#############################
##   triEygenP11: Extension operator for y from Jacobi P^(1,1) basis into T = JacobiTriangle(1,0,0,Triangle([0.0,1.0],[0.0,0.0],[1.0,1.0]))
####
function triEygenP11(N::Integer)
    Ey = BlockBandedMatrix(Zeros(sum(1:N),N), 1:N,Ones{Int}(N), (0,0))
        for n=1:N
                view(Ey, Block(n,n))[1,1] = (-1)^(n+1)
        end
    return Ey
end
#############################
##   triExgenP10: Extension operator for x from Jacobi P^(1,0) basis to canonical Jacobi triangle
####
function triExgenP01(N::Integer)
    Ex = BlockBandedMatrix(Zeros(sum(1:N), N), 1:N, Ones{Int}(N), (0,0))
        for n=1:N
            view(Ex, Block(n,n))[1,1] = 1
        end
    return Ex
end
