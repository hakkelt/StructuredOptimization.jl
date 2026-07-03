using Documenter, StructuredOptimization, 
LinearAlgebra, DSP, FFTW, AbstractOperators, ProximalAlgorithms

makedocs(
  modules = [StructuredOptimization,ProximalAlgorithms],
  format = Documenter.HTML(),
  # Phase 0.3: run every docstring/doc code block as a doctest in CI.
  doctest = true,
  sitename = "StructuredOptimization",
  authors = "Niccolò Antonello and Lorenzo Stella",
  pages = [
  "Home"                  => "index.md",
  "Quick Tutorial Guide"  => "tutorial.md",
  "Expressions"           => "expressions.md",
  "Functions"             => "functions.md",
  "Solvers"               => "solvers.md",
  "Demos"                 => "demos.md",
  ],
)

deploydocs(
  repo   = "github.com/kul-forbes/StructuredOptimization.jl.git",
  target = "build",
)
