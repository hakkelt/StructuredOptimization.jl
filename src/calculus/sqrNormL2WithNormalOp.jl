# squared L2 norm (times a constant, or weighted) precomposed with an operator

"""
    SqrNormL2WithNormalOp(L::LinearOperator, λ=1)

With a nonnegative scalar `λ`, return the squared Euclidean norm
```math
f(x) = \\tfrac{λ}{2}\\|L * x\\|^2.
```
With a nonnegative array `λ`, return the weighted squared Euclidean norm
```math
f(x) = \\tfrac{1}{2}∑_i λ_i y_i^2 where y = L * x.
```

This is a special case of the more general `Precompose(SqrNormL2(), L, 1, 0)` operator,
where `L` is a linear operator, and only the gradient is needed, not the proximal operator.
The gradient of the precomposed squared norm is
```math
\\nabla f(x) = Lᴴ * L * x,
```
and in many cases, there is an optimized implementation of the normal operator `Lᴴ * L`
that makes the compution of the gradient much faster than the naive implementation.

A notable drawback of this method is that gradient! does not return the
squared norm of `L * x`, but rather the squared norm of `Lᴴ * L * x` (i.e. the
squared norm of the gradient). Most algorithms, however, tolerate this
difference, and it is much faster to compute.
"""
struct SqrNormL2WithNormalOp{T, SC, L <: AbstractOperator, L2 <: AbstractOperator}
    A::L
    # Normal operator used for the gradient. For scalar λ it is AᴴA (the weight is
    # applied afterwards); for array λ it is the *weighted* normal operator
    # Aᴴ·diag(λ)·A, so the gradient Aᴴ·diag(λ)·A·x is computed in one mul!.
    AᴴA::L2
    lambda::T
    function SqrNormL2WithNormalOp(A, lambda)
        @assert A isa AbstractOperator
        @assert is_linear(A)
        if any(lambda .< 0)
            error("coefficients in λ must be nonnegative")
        end
        # Strong convexity of x ↦ ½‖diag(√λ)·A·x‖² needs a positive weight *and* an
        # injective operator (full column rank), otherwise the null space of A is flat.
        strongly_convex = all(lambda .> 0) && is_full_column_rank(A)
        if lambda isa AbstractArray
            W = AbstractOperators.DiagOp(AbstractOperators.codomain_type(A), size(A, 1), lambda)
            AᴴA = A' * W * A
        else
            AᴴA = A' * A
        end
        return new{typeof(lambda), strongly_convex, typeof(A), typeof(AᴴA)}(A, AᴴA, lambda)
    end
end

is_convex(::Type{<:SqrNormL2WithNormalOp}) = true
is_smooth(::Type{<:SqrNormL2WithNormalOp}) = true
is_separable(::Type{<:SqrNormL2WithNormalOp}) = true
is_generalized_quadratic(::Type{<:SqrNormL2WithNormalOp}) = true
is_strongly_convex(::Type{<:SqrNormL2WithNormalOp{T, SC}}) where {T, SC} = SC

SqrNormL2WithNormalOp(A) = SqrNormL2WithNormalOp(A, 1)

function (f::SqrNormL2WithNormalOp{S})(x) where {S <: Real}
    y = f.A * x
    return f.lambda / real(eltype(y))(2) * norm(y)^2
end

function (f::SqrNormL2WithNormalOp{<:AbstractArray})(x)
    y = f.A * x
    R = real(eltype(y))
    sqnorm = R(0)
    for k in eachindex(y)
        sqnorm += f.lambda[k] * abs2(y[k])
    end
    return sqnorm / R(2)
end

function gradient!(y, f::SqrNormL2WithNormalOp{<:Real}, x)
    R = real(eltype(y))
    mul!(y, f.AᴴA, x)
    sqnx = R(0)
    for k in eachindex(y)
        y[k] *= f.lambda
        sqnx += abs2(y[k])
    end
    return f.lambda / R(2) * sqnx
end

function gradient!(y, f::SqrNormL2WithNormalOp{<:AbstractArray}, x)
    R = real(eltype(y))
    # f.AᴴA is the weighted normal operator Aᴴ·diag(λ)·A, so this is exactly the
    # gradient ∇f(x) = Aᴴ·diag(λ)·A·x (weights applied in the codomain, not the domain).
    mul!(y, f.AᴴA, x)
    sqnx = R(0)
    for k in eachindex(y)
        sqnx += abs2(y[k])
    end
    return sqnx / R(2)
end
