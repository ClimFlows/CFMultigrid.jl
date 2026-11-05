# CFMultigrid

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ClimFlows.github.io/CFMultigrid.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ClimFlows.github.io/CFMultigrid.jl/dev/)
[![Build Status](https://github.com/ClimFlows/CFMultigrid.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ClimFlows/CFMultigrid.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ClimFlows/CFMultigrid.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ClimFlows/CFMultigrid.jl)

CFMultigrid is a matrix-free geometric multigrid solver for the Poisson and
Helmholtz equations on structured Cartesian grids. It is written in pure Julia
and is used by the [ClimFlows](https://github.com/ClimFlows) atmospheric models.

## Features

- 2-D (single vertical layer) and 3-D Cartesian grids.
- Automatic two-way coarsening of the grid hierarchy.
- `Poisson`, `PoissonNonUniform` (anisotropic coefficients) and `Helmholtz`
  operators.
- `NEUMANN` and `DIRICHLET` boundary conditions.
- Weighted `Jacobi` and vertical `LineRelaxation` smoothers.
- V-cycles driven by a small `Param` object.

## Installation

CFMultigrid is distributed through the ClimFlows registry:

```julia
using Pkg
Pkg.Registry.add(url = "https://github.com/ClimFlows/JuliaRegistry.git")
Pkg.add("CFMultigrid")
```

For development, clone the repository and use `Pkg.develop`:

```julia
using Pkg
Pkg.develop(path = "/path/to/CFMultigrid.jl")
```

## Quick start

```julia
using CFMultigrid

# Number of grid cells in x, y and z, plus the number of halo cells.
shape = (64, 64, 64)
nhalo = 1

# Pre-smoothing steps, post-smoothing steps, coarsest-level smoothing steps,
# maximum number of V-cycles, tolerance and smoother.
param = Param(1, 2, 20, 20, 1e-9, Jacobi(0.85))

# Build the multigrid hierarchy for the Poisson operator with Neumann BCs.
mg = setup(shape, nhalo, Poisson, NEUMANN)

# Fill the right-hand side of the finest grid with a dipole.
b = mg[1].b
b[1 + div(64, 3), 1 + div(64, 3), 1] =  1
b[64 - div(64, 3), 64 - div(64, 3), 1] = -1

# Solve; the initial guess is mg[1].x (zero by default).
res, nite = solve!(mg, param, VCYCLE; verbose = false)
println("res = ", res, " after ", nite, " V-cycles")
```

This produces:

```text
res = 7.647626690763077e-10 after 10 V-cycles
```

After the solve, the solution is stored in `mg[1].x`. The `mg` object is a
tuple of grids ordered from the finest to the coarsest level; `show(mg)` gives
an overview of the hierarchy.

## Grids and operators

`setup` builds the hierarchy from the finest grid down to a `2 x 2 x 2` (or
`2 x 2`) coarsest grid. Grid dimensions are expressed in cells and are
unit-spaced; use power-of-two sizes for a clean coarsening. A value of
`nz = 1` selects the 2-D (`XY`) variant, while `nz > 1` selects the 3-D (`XYZ`)
variant.

The available operators are:

| Operator | Description |
| --- | --- |
| `Poisson` | Constant-coefficient Poisson operator (7-point in 3-D, 5-point in 2-D). |
| `PoissonNonUniform` | Poisson operator with anisotropic coefficients `cxx`, `cyy`, `czz`. |
| `Helmholtz` | Poisson operator with an additional `-coef * x` term, i.e. `(∇² - coef) x = b`. |

Keyword arguments accepted by `setup` are forwarded to the operator:

```julia
# Anisotropic Poisson with line relaxation in the vertical direction.
param = Param(1, 2, 20, 20, 1e-9, LineRelaxation(8, Float64))
mg = setup((128, 128, 8), 1, PoissonNonUniform, NEUMANN;
           cxx = 1.0, cyy = 1.0, czz = 100.0)

# Helmholtz operator: (∇² - coef) x = b, using unit coefficients in x and y.
mg = setup((64, 64, 1), 1, Helmholtz, NEUMANN;
           cxx = 1.0, cyy = 1.0, czz = 100.0, coef = 1.0)
```

The operators accept an optional `msk` keyword argument (`Int8` array, with
`1` for active cells and `0` for excluded cells). It is not yet usable through
`setup`, because the same mask is currently forwarded to every grid level; see
[issue #17](https://github.com/ClimFlows/CFMultigrid.jl/issues/17).

## API overview

- `setup(shape, nhalo, operator, bc; kwargs...)` builds the multigrid hierarchy
  and returns a tuple of grids, finest first.
- `Param(npre, npost, ndeepest, maxite, tol, smoother)` configures the cycle:
  number of pre- and post-smoothing steps, number of smoothing steps on the
  coarsest level, maximum number of V-cycles, relative L2-residual tolerance
  and smoother.
- `solve!(mg, param, VCYCLE; verbose = true)` runs V-cycles until the relative
  L2 residual `‖r‖₂ / ‖b‖₂` drops below `param.tol` or `param.maxite` is
  reached. It returns `(res, nite)`, where `res` is that relative L2 residual.
  For example, `tol = 1e-9` requests `‖r‖₂ / ‖b‖₂ <= 1e-9`.
- Smoothers: `Jacobi(omega)` and `LineRelaxation(n, T)`.
- Boundary conditions: `NEUMANN` and `DIRICHLET`.
- Low-level transfer operators `restriction` and `prolongation` are exported
  for advanced use.

## Testing

```julia
using Pkg
Pkg.test("CFMultigrid")
```

The test suite exercises the Poisson, non-uniform Poisson and Helmholtz
solvers in 2-D and 3-D.

## Documentation

The documentation is available at
[ClimFlows.github.io/CFMultigrid.jl](https://ClimFlows.github.io/CFMultigrid.jl/stable/).
It is built from the docstrings in `src/`.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the release history.

## License

CFMultigrid is released under the MIT license; see [LICENSE](LICENSE).

## Contributing

Contributions are welcome through pull requests. Changes merged to `main` are
tested on Julia 1.10, 1.11 and 1.12. When the `version` field of
`Project.toml` changes on `main`, CI automatically tags the corresponding
release and registers it in the ClimFlows registry.
