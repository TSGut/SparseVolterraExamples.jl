module SparseVolterraExamples

    using ApproxFun, NLsolve, BandedMatrices, BlockBandedMatrices, LinearAlgebra, BlockArrays, MultivariateOrthogonalPolynomials, Test, Plots

    include("triIntegroDiff.jl")
    include("triQEgen.jl")
    include("triVolterraFullKernel.jl")
    include("triVolterraMonKernel.jl")

    export triNonLinearVolterraObjective,
        triNonLinearIntegroDiffVolterraObjective,
        pochhammer,
        alternI,
        EvalJPat1mat,
        DirectEvalLHSP10at0,
        DirectEvalLHSP10atPrime0,
        DirectEvalLHSP10atPPrime0,
        DirectEvalLHSP10atPPPrime0,
        DirectEvalLHSP10atPPPPrime0,
        triQgen,
        triEygenP01,
        triQEygenP01,
        WLoweringP01P00,
        reflectPabtoPba,
        WLoweringP01Pm11,
        ReflectorP01,
        triEygenP00,
        triEygenP11,
        triExgenP01,
        triVolterraEQ1FullKernelSolver,
        triVolterraEQ2FullKernelSolver,
        triVolterraFullKernelFunP01,
        triVolterraFullKernelOpP01,
        triVolterraFullKernelOpP01,
        triBandedSgen,
        triVolterraEQ1Solver,
        triVolterraEQ2Solver,
        triVolterraFunP01,
        triVolterraOpP01,
        triJQEJbasislexP01,
        triJQEJbasisgridgenP01,
        triNoKernelDirect

end # module
