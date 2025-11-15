struct Expression{N,A<:AbstractOperator} <: AbstractExpression
  x::NTuple{N,Variable}
  L::A
  function Expression(x::NTuple{N,Variable}, L::A) where {N,A<:AbstractOperator}
    # checks on L
    ndoms(L,1) > 1 && throw(ArgumentError(
      "Cannot create expression with LinearOperator with `ndoms(L,1) > 1`"
     ))
    #checks on x
    szL = size(L,2)
    szx = size.(x)
    check_sz = length(szx) == 1 ? szx[1] != szL : szx != szL
    check_sz && throw(ArgumentError(
      "Size of the operator domain $(size(L, 2)) must match size of the variable $(size.(x))"
     ))
    dmL = domain_type(L)
    dmx = eltype.(x)
    check_dm = length(dmx) == 1 ? dmx[1] != dmL : dmx != dmL
    check_dm && throw(ArgumentError(
      "Type of the operator domain $(domain_type(L)) must match type of the variable $(eltype.(x))"
     ))
    new{N,A}(x,L)
  end
end

struct AdjointExpression{E <: AbstractExpression} <: AbstractExpression
  ex::E
end

import Base: adjoint, show

adjoint(ex::AbstractExpression) = AdjointExpression(convert(Expression,ex))
adjoint(ex::AdjointExpression) = ex.ex

function show(io::IO, ex::Expression)
  if length(ex.x) == 1
    print(io, AbstractOperators.fun_name(ex.L), " * ", ex.x[1])
  else
    print(io, AbstractOperators.fun_name(ex.L), " * (", join(ex.x, ", "), ")")
  end
end

include("utils.jl")
include("multiplication.jl")
include("addition.jl")
include("addition_tricky_part.jl")
include("abstractOperator_bind.jl")
