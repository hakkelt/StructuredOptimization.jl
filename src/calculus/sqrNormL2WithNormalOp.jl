# squared L2 norm (times a constant, or weighted) precomposed with an operator

"""
    SqrNormL2WithNormalOp(L::AbstractOperator, λ=1)

With a nonnegative scalar `λ`, return the squared Euclidean norm
```math
f(x) = \\tfrac{λ}{2σ}\\|L * x\\|^2,
```
where `σ` is the adjoint scaling of `L` described below (`σ = 1`, and the factor disappears,
whenever `L'` is the true adjoint of `L`).
With a nonnegative array `λ`, return the weighted squared Euclidean norm
```math
f(x) = \\tfrac{1}{2σ}∑_i λ_i y_i^2 where y = L * x.
```

This is a special case of the more general `Precompose(SqrNormL2(), L, 1, 0)` operator,
where `L` is a linear operator, and only the gradient is needed, not the proximal operator.
The gradient of the precomposed squared norm is
```math
\\nabla f(x) = Lᴴ * L * x,
```
and in many cases, there is an optimized implementation of the normal operator `Lᴴ * L`
that makes the computation of the gradient much faster than the naive implementation.

`L` may be affine (an `AffineAdd`, as produced by `ls(A*x - b)`): writing `L*x = A*x + d`,
the normal operator carries the displacement `Aᴴd` automatically (`Lᴴ*L*x = AᴴA*x + Aᴴd`
when `L*0 = d`), so `gradient!` computes the correct gradient in a single pass.

`gradient!` returns the function value `f(x)`, as `ProximalCore.value_and_gradient!`
requires, recovered from the gradient without a second application of `L`.

# Adjoint scaling

`L'` is not always the true adjoint of `L`. A `BACKWARD`-normalized DFT, for instance, has
`L' = L⁻¹ = Lᴴ/N`: the pair is off by a positive scalar `σ` defined by
```math
\\mathrm{Re}⟨L u, L u⟩ = σ \\, \\mathrm{Re}⟨u, (L'L) u⟩ .
```
Since `Lᴴ*L*x` (as actually computed from `L'*L`) is then `1/σ` times the true gradient of
`f`, the value returned alongside it must be scaled the same way for the two to be
consistent — otherwise anything that reads both (a backtracking line search, a printed
objective) is meaningless. `σ` is measured once, at construction, with a single probe
through `L` and `L'L`; with a genuine adjoint it is `1` and every formula above reduces to
the usual one.
"""
struct SqrNormL2WithNormalOp{T, SC, L <: AbstractOperator, L2 <: AbstractOperator, D, R <: Real}
    A::L
    # Normal operator used for the gradient. For scalar λ it is AᴴA (the weight is
    # applied afterwards); for array λ it is the *weighted* normal operator
    # Aᴴ·diag(λ)·A, so the gradient Aᴴ·diag(λ)·A·x is computed in one mul!.
    AᴴA::L2
    lambda::T
    # `Aᴴd`: the normal operator's displacement (`AᴴA * 0`), taken through the same,
    # possibly weighted, operator `gradient!` uses, or `nothing` when `A` is purely
    # linear (the overwhelmingly common case), so the per-gradient correction is
    # skipped entirely rather than paying a dot product with zeros.
    Aᴴd::D
    # The constant term of the quadratic, `‖d‖²/(2σ)` (weighted by λ when λ is an array).
    half_sqnorm_d::R
    # `1/σ`, the adjoint scaling of `A` (see the docstring); `1` for a true adjoint pair.
    inv_scaling::R
    function SqrNormL2WithNormalOp(A, lambda)
        @assert A isa AbstractOperator
        @assert is_linear(A)
        if any(lambda .< 0)
            error("coefficients in λ must be nonnegative")
        end
        # Strong convexity of x ↦ ½‖diag(√λ)·A·x‖² needs a positive weight *and* an
        # injective operator (full column rank), otherwise the null space of A is flat.
        strongly_convex = all(lambda .> 0) && is_full_column_rank(A)
        # Built unweighted, purely to measure the adjoint scaling below: that scaling is a
        # property of the (A, A') pair alone and is unaffected by inserting a Hermitian,
        # positive weight between them.
        pureAᴴA = A' * A
        if lambda isa AbstractArray
            W = AbstractOperators.DiagOp(AbstractOperators.codomain_type(A), size(A, 1), lambda)
            AᴴA = A' * W * A
        else
            AᴴA = lambda == 1 ? pureAᴴA : lambda * pureAᴴA
        end
        # `A * 0` is the displacement `d` of an affine `A` (zero for a purely linear one);
        # `AᴴA * 0` is `Aᴴd` taken through the very operator `gradient!` uses, so the
        # constants cannot drift from it.
        z = AbstractOperators.allocate_in_domain(A)
        fill!(z, 0)
        d = A * z
        has_displacement = !iszero(d)
        Aᴴd = has_displacement ? AᴴA * z : nothing
        inv_scaling = _inv_adjoint_scaling(A, pureAᴴA, z, d, has_displacement ? pureAᴴA * z : nothing)
        R_ = typeof(inv_scaling)
        half_sqnorm_d = has_displacement ? R_(_weighted_sqnorm(lambda, d) * inv_scaling / 2) : zero(R_)
        return new{typeof(lambda), strongly_convex, typeof(A), typeof(AᴴA), typeof(Aᴴd), R_}(
            A, AᴴA, lambda, Aᴴd, half_sqnorm_d, inv_scaling
        )
    end
end

_weighted_sqnorm(lambda::Real, d) = lambda * real(dot(d, d))
function _weighted_sqnorm(lambda::AbstractArray, d)
    R = real(eltype(d))
    sqnorm = R(0)
    for k in eachindex(d)
        sqnorm += lambda[k] * abs2(d[k])
    end
    return sqnorm
end

# `σ` from the docstring, as `1/σ`: `Re⟨A u, A u⟩ / Re⟨u, (A'A) u⟩` for a probe `u`, with the
# displacement of an affine `A` subtracted so that only the linear parts are compared.
#
# The probe is the constant vector, which is deterministic (no RNG dependency, so the value a
# solver prints does not move between runs) and is annihilated by no operator this is used
# with. Should it nevertheless land in the null space, `Aᴴd` — nonzero exactly when there is a
# displacement to correct — is tried next; if that fails too the scaling is left at 1, which is
# the behaviour of a true adjoint pair.
function _inv_adjoint_scaling(A, AᴴA, z, d, Aᴴd)
    R = real(eltype(z))
    u = similar(z)
    for probe in 1:2
        if probe == 1
            fill!(u, one(eltype(z)))
        elseif Aᴴd !== nothing
            copyto!(u, Aᴴd)
        else
            break
        end
        Au = A * u
        w = AᴴA * u
        # strip the affine displacement: `A u = A_lin u + d` and `(A'A) u = (A'A)_lin u + Aᴴd`
        if Aᴴd !== nothing
            Au = Au .- d
            w = w .- Aᴴd
        end
        num = real(dot(Au, Au))
        den = real(dot(u, w))
        isfinite(num) && isfinite(den) && den > 0 && return R(den / num)
    end
    return one(R)
end

is_convex(::Type{<:SqrNormL2WithNormalOp}) = true
is_smooth(::Type{<:SqrNormL2WithNormalOp}) = true
is_separable(::Type{<:SqrNormL2WithNormalOp}) = true
is_generalized_quadratic(::Type{<:SqrNormL2WithNormalOp}) = true
is_strongly_convex(::Type{<:SqrNormL2WithNormalOp{T, SC}}) where {T, SC} = SC

SqrNormL2WithNormalOp(A) = SqrNormL2WithNormalOp(A, 1)

function (f::SqrNormL2WithNormalOp)(x)
    y = f.A * x
    return _weighted_sqnorm(f.lambda, y) * f.inv_scaling / 2
end

function gradient!(y, f::SqrNormL2WithNormalOp, x)
    mul!(y, f.AᴴA, x)
    v = real(dot(x, y)) / 2
    if f.Aᴴd !== nothing
        v += real(dot(x, f.Aᴴd)) / 2 + f.half_sqnorm_d
    end
    return v
end
