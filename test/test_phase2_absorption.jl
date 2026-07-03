# Phase 2.1 — property tests for the canonical absorption transform
# `merge_function_with_operator(op, f, disp, λ)`, which must satisfy
#     absorbed_f(x) ≈ λ · f(op * x + disp)
# for every absorption case (eye / diagonal / AAᴴ-diagonal / general linear /
# nonlinear). This is the value-level invariant that the Phase 1 displacement/λ
# bugs (1.1–1.4) all violated.

const SO2 = StructuredOptimization
const merge_fo = StructuredOptimization.merge_function_with_operator

@testset "Phase 2.1 absorption property" begin
    Random.seed!(200)

    # helper: absorbed(x) ≈ λ * f(op*x + disp)
    function check_absorption(op, f, disp, λ; cplx=false)
        g = merge_fo(op, f, disp, λ)
        for _ in 1:3
            x = cplx ? randn(ComplexF64, size(op, 2)) : randn(size(op, 2))
            expected = λ * f(op * x .+ disp)
            @test abs(g(x) - expected) < 1e-9 * (1 + abs(expected))
        end
    end

    # eye operator with displacement (non-SqrNormL2 function)
    @testset "eye" begin
        op = AbstractOperators.Eye(Float64, (4,))
        check_absorption(op, NormL1(), randn(4), 2.5)
        check_absorption(op, NormL1(), zeros(4), 1.0)
    end

    # diagonal operator: NormL1 keeps displacement; SqrNormL2 folds into the weight
    # (valid only at zero displacement, which is how it is reached in prepare).
    @testset "diagonal" begin
        op = DiagOp([2.0, -3.0, 4.0, 0.5])
        check_absorption(op, NormL1(), randn(4), 1.7)
        check_absorption(op, SqrNormL2(1.5), zeros(4), 2.0)
    end

    # AAᴴ-diagonal operator (DFT): AAᴴ = N·I, hit via the Precompose branch
    @testset "AAc-diagonal" begin
        xv = Variable(8)
        op = SO2.operator(fft(xv))          # DFT: ℝ^8 -> ℂ^8, AAᴴ = 8·I
        check_absorption(op, SqrNormL2(), randn(ComplexF64, 8), 1.3)
    end

    # general (non-square, non-AAᴴ-diagonal) linear operator
    @testset "general linear" begin
        op = MatrixOp(randn(6, 4))
        check_absorption(op, NormL1(), randn(6), 0.9)
        check_absorption(op, SqrNormL2(), randn(6), 2.2)
    end

    # nonlinear operator (sin): PrecomposeNonlinear with AffineAdd displacement
    @testset "nonlinear" begin
        xv = Variable(5)
        op = SO2.operator(sin(xv))
        check_absorption(op, SqrNormL2(), randn(5), 3.1)
    end
end
