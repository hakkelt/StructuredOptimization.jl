module StructuredOptimization

using LinearAlgebra
using RecursiveArrayTools
using ProximalCore
using AbstractOperators
using ProximalOperators
using ProximalAlgorithms
using Combinatorics: permutations, powerset
using OperatorCore

ProximalAlgorithms.value_and_gradient(f, x) = begin
  y, fy = gradient(f, x)
  return fy, y
end
abstract type AbstractExpression end

include("syntax/variable.jl")
include("syntax/expressions/expression.jl")
include("syntax/terms/term.jl")

const TermOrExpr =  Union{Term,AbstractExpression}

include("calculus/precomposeNonlinear.jl") # TODO move to ProximalOperators?
include("calculus/sqrNormL2WithNormalOp.jl")

# problem parsing
include("solvers/terms_extract.jl")
include("solvers/terms_properties.jl")
include("solvers/parse.jl")

# solver calls
include("solvers/build_solve.jl")
include("solvers/minimize.jl")

end
