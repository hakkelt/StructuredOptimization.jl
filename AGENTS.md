# AGENTS.md — StructuredOptimization.jl

StructuredOptimization.jl is a high-level Julia interface for composite optimization problems of the form

    minimize  f(Ax) + g(x)

It provides an algebraic syntax for building expressions and problems from `Variable`s, then dispatches to first-order algorithms from **ProximalAlgorithms.jl**.

## Architecture

```
Variable → AbstractExpression → Term → problem() → solve()
```

| Layer | Files | Role |
|---|---|---|
| Syntax | `src/syntax/variable.jl`, `src/syntax/expressions/`, `src/syntax/terms/` | Build operator graphs |
| Calculus | `src/calculus/` | Custom proximal operators: `SqrNormL2WithNormalOp`, `precomposeNonlinear` |
| Solvers | `src/solvers/` | Extract terms, parse problem structure, dispatch algorithms |

Key solver files:
- `terms_extract.jl` — extract variables, operators, affines, functions from a `Term`
- `terms_properties.jl` — classify terms (proximable, smooth, etc.)
- `parse.jl` — match problem structure to algorithm assumptions
- `build_solve.jl` — `solve()`, `print_diagnostics()`, `suggest_algorithm()`
- `minimize.jl` — `@minimize` macro

**Dependencies**: `AbstractOperators.jl`, `ProximalOperators.jl`, `ProximalAlgorithms.jl`, `ProximalCore.jl` are dev'd locally via `test/Project.toml` `[sources]`, pointing at sibling checkouts (`../../AbstractOperators`, `../../ProximalAlgorithms.jl`, etc.). Those checkouts may be on feature branches — check `git -C <path> branch --show-current` rather than assuming a branch name, since it changes over time.

Two declared deps are not referenced in `src/` today and are kept intentionally:
- `DifferentiationInterface` — reserved for Phase 5 (differentiable solvers / unrolling); `ProximalAlgorithms` already routes autodiff through it.
- `AbstractFFTs` — the generic `fft`/`rfft` interface that `FFTWOperators` and the `import FFTW: fft` bindings build on; retained so the FFT expression bindings resolve against a stable interface package.
`DSP`/`FFTW` are used (function-name imports in `syntax/expressions/abstractOperator_bind.jl`).

## Testing Conventions

- Test files are standalone modules (prefix `test_`) included from `test/runtests.jl`
- Deterministic tests: `Random.seed!(0)` is set globally in `test/runtests.jl`; individual test files may reset with `Random.seed!(n)` for isolated seeds
- Prefer **PANOCplus** for optimization tests — PANOC and ZeroFPR hit "stepsize too small" on many problems and produce unreliable results. Only test PANOC/ZeroFPR when testing solver dispatch, and mark known-failing convergence checks as `@test_broken`
- `Aqua.jl` runs in `runtests.jl` with `ambiguities=false, piracies=false, persistent_tasks=false` at the top level, plus separate `broken=true` checks for ambiguities/persistent_tasks and an explicit piracy allowlist (`treat_as_own`) for the ProximalAlgorithms/ProximalOperators methods this package legitimately extends

### Algorithm Selection Guide
| Problem type | Recommended solver |
|---|---|
| `f(Ax) + g(x)`, f smooth | `PANOCplus` |
| Pure proximal (`g(x)` only) | `FastForwardBackward` |
| Comparison across solvers | use `PANOCplus` and `FastForwardBackward`; add `ZeroFPR` only if testing dispatch |
| Avoid for convergence tests | `PANOC` — unreliable stepsize; `ZeroFPR` — sometimes hits stepsize-too-small |

## Development Workflow

### Environment Setup
The test environment is separate from the package environment:
```sh
cd test/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Running Tests

**Full test suite**:
```sh
julia --project=test -e '
  using StructuredOptimization, AbstractOperators, DSPOperators, FFTWOperators
  using ProximalOperators, ProximalAlgorithms, RecursiveArrayTools
  using LinearAlgebra, Random, DSP, FFTW, Test
  include("test/runtests.jl")
'
```

**Single test file** (from the package root):
```sh
julia --project=test -e '
  using StructuredOptimization, AbstractOperators, ProximalOperators, ProximalAlgorithms
  using RecursiveArrayTools, LinearAlgebra, Random, Test
  Random.seed!(0)
  include("test/test_usage.jl")
'
```

### Coverage
Use `LocalCoverage.jl` to measure and view test coverage:
```sh
julia --project=test -e '
  using LocalCoverage
  cov = generate_coverage("StructuredOptimization")
  LocalCoverage.html_coverage(cov; dir="coverage_html")
'
```
`generate_coverage` runs the test suite with `--code-coverage=user` and drops `*.jl.<pid>.cov` files next to each source file — remove them (`find . -name '*.cov' -delete`) once done, they are generated artifacts and should not be committed.

### Formatting
- This project uses **Runic.jl** for formatting
- Install: `julia --project=@runic --startup-file=no -e 'using Pkg; Pkg.add("Runic")'`
- Format: `julia --project=@runic --startup-file=no -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace src/`
- Format before committing

## Known Issues / Broken Tests

| Test | Status | Root cause |
|---|---|---|
| PANOC lasso/box/NNLS convergence in `test_usage.jl` | `@test_broken` | Upstream PANOC stepsize-too-small bug in ProximalAlgorithms.jl |
| Aqua ambiguities (`Base.:+`, `Base.:<=`, `Base.:>=`) | `@test_broken`, excluded | Ambiguities from this package's operator overloads |

## Package Structure

```
src/
  StructuredOptimization.jl   # module entry; SeparableSum ArrayPartition dispatch
  syntax/
    variable.jl               # Variable type, ~x dereference, get_name
    expressions/
      expression.jl           # AbstractExpression, operator(), affine(), variables()
      addition.jl             # Usum_op, expression + expression
      addition_tricky_part.jl # add_missing_vars, multi-variable sum support
      ...
    terms/
      term.jl                 # Term type, ls(), norm(), smooth(), ...
  calculus/
    precomposeNonlinear.jl
    sqrNormL2WithNormalOp.jl  # normalop_ls
  solvers/
    terms_extract.jl
    terms_properties.jl
    parse.jl
    build_solve.jl
    minimize.jl
test/
  runtests.jl
  test_variables.jl
  test_expressions.jl
  test_AbstractOp_binding.jl
  test_terms.jl
  test_proxstuff.jl
  test_problem.jl
  test_build_minimize.jl
  test_usage_small.jl
  test_usage.jl
```
