struct Term{T1 <: Real, T2, T3 <: AbstractExpression}
    lambda::T1
    f::T2
    A::T3
    repr::Union{String,Nothing}
end

function Term(lambda, f, ex::AbstractExpression)
    return Term(lambda,f,ex,nothing)
end

function Term(f, ex::AbstractExpression)
    A = convert(Expression,ex)
    Term(one(real(codomainType(affine(A)))),f, A)
end

function Term(f, ex::AbstractExpression, repr::String)
    A = convert(Expression,ex)
    Term(one(real(codomainType(affine(A)))),f, A, repr)
end

function Term(t::Term, repr::String)
    Term(t.lambda, t.f, t.A, repr)
end

import Base: ==, show

# Ignore the repr when comparing terms
==(t1::Term, t2::Term) = t1.lambda == t2.lambda && t1.f == t2.f && t1.A == t2.A

function show(io::IO, t::Term)
    if t.repr !== nothing
        print(io, t.repr)
    else
        print(io, t.lambda, " * ", t.f, "(", t.A, ")")
    end
end

function show(io::IO, t::NTuple{N,Term}) where {N}
    for i in 1:N
        show(io, t[i])
        if i < N
            print(io, " + ")
        end
    end
end

# Operations

# Define sum of terms simply as their vcat

import Base: +

(+)(a::Term,b::Term) = (a,b)
(+)(a::NTuple{N,Term},b::Term)    where {N} = (a...,b)
(+)(a::Term,b::NTuple{N,Term})    where {N} = (a,b...)
(+)(a::NTuple{N,Term},::Tuple{}) where {N} = a
(+)(::Tuple{},b::NTuple{N,Term}) where {N} = b
(+)(a::NTuple{N,Term},b::NTuple{M,Term}) where {N,M} = (a...,b...)

# Define multiplication by constant

import Base: *

function (*)(a::T1, t::Term{T,T2,T3}) where {T1<:Real, T, T2, T3}
    coeff = *(promote(a,t.lambda)...)
    Term(coeff, t.f, t.A)
end

function (*)(a::T1, t::T2) where {T1<:Real, N, T2 <: Tuple{Vararg{<:Term,N}} }
    return a.*t
end

# Properties

variables(t::Term) = variables(t.A)
operator(t::Term) = operator(t.A)
affine(t::Term) = affine(t.A)
displacement(t::Term) = displacement(t.A)

#importing properties from ProximalOperators
import ProximalCore:
              is_affine_indicator,
              is_cone_indicator,
              is_convex,
              is_generalized_quadratic,
              is_proximable,
              is_quadratic,
              is_separable,
              is_set_indicator,
              is_singleton_indicator,
              is_smooth,
              is_locally_smooth,
              is_strongly_convex

is_func_f = [
        :is_set_indicator,
        :is_singleton_indicator,
        :is_smooth,
        :is_locally_smooth,
       ]

for f in is_func_f
    @eval begin
        import ProximalCore: $f
        $f(t::Term) = $f(t.f)
        $f(t::NTuple{N,Term}) where {N} = all($f.(t))
    end
end

#importing properties from AbstractOperators
is_op_f = [:is_linear,
        :is_eye,
        :is_null,
        :is_diagonal,
        :is_AcA_diagonal,
        :is_AAc_diagonal,
        :is_orthogonal,
        :is_invertible,
        :is_full_row_rank,
        :is_full_column_rank,
        :is_sliced
       ]

for f in is_op_f
    @eval begin
        import AbstractOperators: $f
        $f(t::Term) = $f(operator(t))
        $f(t::NTuple{N,Term}) where {N} = all($f.(t))
    end
end

is_affine_indicator(t::Term) = is_affine_indicator(t.f) && is_linear(t)
is_cone_indicator(t::Term) = is_cone_indicator(t.f) && is_linear(t)
is_convex(t::Term)    = is_convex(t.f) && is_linear(t)
is_quadratic(t::Term) = is_quadratic(t.f) && is_linear(t)
is_generalized_quadratic(t::Term) = is_generalized_quadratic(t.f) && is_linear(t)
is_strongly_convex(t::Term) = is_strongly_convex(t.f) && is_full_column_rank(operator(t.A))
is_separable(t::Term) = is_separable(t.f) && is_diagonal(operator(t.A))

include("proximalOperators_bind.jl")

# other stuff, to make Term work with iterators
import Base: iterate, isempty
iterate(t::Term, state = true) = state ? (t, false) : nothing
isempty(t::Term) =  false
