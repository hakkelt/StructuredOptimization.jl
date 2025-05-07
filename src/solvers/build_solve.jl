const ForwardBackwardSolver = ProximalAlgorithms.IterativeAlgorithm

"""
	parse_problem(terms::Tuple, solver::ForwardBackwardSolver)

Takes as input a tuple containing the terms defining the problem and the solver.

Returns a tuple containing the optimization variables and the problem terms
to be fed into the solver.

# Example

```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem( ls(A*x - b ) , norm(x) <= 1 );

julia> StructuredOptimization.parse_problem(p, PANOCplus());
```
"""
function parse_problem(terms::NTuple{N,StructuredOptimization.Term}, algorithm::T, return_partial::Bool = false) where {N,T <: ForwardBackwardSolver}
    assumptions = ProximalAlgorithms.get_assumptions(algorithm)
    variables = StructuredOptimization.extract_variables(terms)
    remaining_terms = terms
    kwargs = Dict{Symbol, Any}()
    for assumption in assumptions
        for term_selection in reverse(collect(powerset(remaining_terms, 1)))
            term_selection = tuple(term_selection...)
            preparation_result = StructuredOptimization.prepare(term_selection, assumption, variables)
            if preparation_result !== nothing
                term_selection = collect(term_selection)
                remaining_terms = setdiff(remaining_terms, term_selection)
                push!(kwargs, preparation_result...)
                break
            end
        end
        if isempty(remaining_terms)
            return algorithm, kwargs, variables
        end
    end
    return return_partial ? (kwargs, remaining_terms) : nothing
end

function print_diagnostics(terms::NTuple{N,StructuredOptimization.Term}, algorithm::T) where {N,T <: ForwardBackwardSolver}
    kwargs, remaining_terms = parse_problem(terms, algorithm, true)
    print("The algorithm $algorithm assumes problem of form: ")
    show(ProximalAlgorithms.get_assumptions(algorithm))
    if !isempty(kwargs)
        println("Successfully prepared the following terms:")
        for (key, value) in kwargs
            println(" - $key: $value")
        end
    end
    println("The following terms could not be prepared:")
    for term in remaining_terms
        println(" - $term")
    end
end

function parse_problem(terms::NTuple{N,StructuredOptimization.Term}) where {N}
    for algorithm in ProximalAlgorithms.get_algorithms()
        result = parse_problem(terms, algorithm)
        if result !== nothing
            return result
        end
    end
    return nothing
end

function suggest_algorithm(terms::NTuple{N,StructuredOptimization.Term}) where {N}
    suitable_algs = []
    for algorithm in ProximalAlgorithms.get_algorithms()
        result = parse_problem(terms, algorithm)
        if result !== nothing
            push!(suitable_algs, algorithm)
        end
    end
    return suitable_algs
end

function print_diagnostics(terms::NTuple{N,StructuredOptimization.Term}) where {N}
    best_algorithm, best_algorithm_remaining_terms = nothing, Inf
    for algorithm in ProximalAlgorithms.get_algorithms()
        _, remaining_terms = parse_problem(terms, algorithm, true)
        if length(remaining_terms) < best_algorithm_remaining_terms
            best_algorithm_remaining_terms = length(remaining_terms)
            best_algorithm = algorithm
        end
    end
	println("The closest algorithm to the problem is $best_algorithm")
    print_diagnostics(terms, best_algorithm)
end

export solve

"""
	solve(terms::Tuple, solver::ForwardBackwardSolver)

Takes as input a tuple containing the terms defining the problem and the solver options.

Solves the problem returning a tuple containing the iterations taken and the build solver.

# Example

```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem(ls(A*x - b ), norm(x) <= 1);

julia> solve(p, PANOCplus());

julia> ~x
```
"""
function solve(terms::Tuple, solver::ForwardBackwardSolver)
	result = parse_problem(terms, solver)
	if result === nothing
		print_diagnostics(terms, solver)
		error("Sorry, I cannot parse this problem for solver of type $(solver)")
	end
	_, kwargs, x = result
    x_star, it = solver(; x0 = ~x, kwargs...)
	~x .= x_star isa Tuple ? x_star[1] : x_star
	return x, it
end

function solve(terms::Tuple)
	result = parse_problem(terms)
	if result === nothing
		print_diagnostics(terms)
		error("Sorry, I cannot find a suitable solver for this problem")
	end
	solver, kwargs, x = result
	@show solver
    x_star, it = solver(; x0 = ~x, kwargs...)
	~x .= x_star
	return x, it
end
