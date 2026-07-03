using Documenter, StructuredOptimization,
LinearAlgebra, DSP, FFTW, AbstractOperators, ProximalAlgorithms

makedocs(
  # Only this package's exported symbols are coverage-checked; ProximalAlgorithms
  # docstrings are still rendered where referenced, but we don't require documenting
  # its entire internal API here.
  modules = [StructuredOptimization],
  checkdocs = :exports,
  format = Documenter.HTML(),
  # Phase 0.3: run every docstring/doc code block as a doctest in CI.
  doctest = true,
  sitename = "StructuredOptimization",
  authors = "Niccolò Antonello and Lorenzo Stella",
  pages = [
  "Home"                  => "index.md",
  "Quick Tutorial Guide"  => "tutorial.md",
  "Theory" => [
    "How problems are parsed" => "theory/parsing.md",
  ],
  "Expressions"           => "expressions.md",
  "Functions"             => "functions.md",
  "Solvers"               => "solvers.md",
  "FAQ / Troubleshooting" => "faq.md",
  "Demos"                 => "demos.md",
  ],
)

deploydocs(
  repo   = "github.com/hakkelt/StructuredOptimization.jl.git",
  target = "build",
)
