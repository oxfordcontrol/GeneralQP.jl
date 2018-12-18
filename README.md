# GeneralQP.jl
This package solves is a Julia implementation of [the following paper](https://link.springer.com/article/10.1007/BF01588976)
```
Gill, Philip E., and Walter Murray.
Numerically stable methods for quadratic programming.
Mathematical programming 14.1 (1978): 349-372.
```
i.e. an inertia-controlling active set solver for general (definite/indefinite) *dense* quadratic programs
```
minimize    ½x'Px + q'x
subject to  Ax ≤ b
```
given an initial feasible point `x`. 

To avoid further restrictions regarding the initial point, an artificial constraints approach is taken as described in [QPOPT's 1.0 User manual, Section 3.2](https://web.stanford.edu/group/SOL/guides/qpopt.pdf).

## Installation
The solver can be installed by running
```
Pkg.add("https://github.com/oxfordcontrol/TRS.jl")
```
## Usage
The solver can be used by calling the function
```
solve(P, q, A, b, x_init; kwargs) -> x
```
with inputs:

* `P::Matrix{T}` is the quadratic of the cost;
* `q::Vector{T}` is the linear cost;
* `A::Matrix{T}` and b::AbstractVector{T} define the constraints; and
* `x_init` is the initial point

keywords (optional):
* `verbosity::Int=1` the verbosity of the solver ranging from `0` (no output) to `2` (most verbose).
* `printing_interval::Int=50`.

and output `x::Vector{T}`, the calculated optimizer.

## Updatable factorizations
This package includes [`UpdatableQR`](https://github.com/oxfordcontrol/GeneralQP.jl/blob/master/src/linear_algebra.jl), an updatable `QR` factorization `F` of a "thin" `n x m` matrix
```
X = F.Q*F.R
```
that allows efficient `O(n^2)` update of the factors when adding/removing rows in the matrix X.

These can be simply performed using
```
add_column!(F::UpdatableQR{T}, a::AbstractVector{T})
add_column!(F::UpdatableQR{T}, idx::Int).
```

Similarly, [`UpdatableHessianLDL`](https://github.com/oxfordcontrol/GeneralQP.jl/blob/master/src/linear_algebra.jl) provides an updatable LDLt factorization F for the projection of the hessian `P` on the range of the working constraints.

The `struct` `UpdatableHessianLDL` is based on `UpdatableQR` and implements functionality for artificial constraints ([QPOPT 1.0 Manual, Section 3.2](https://web.stanford.edu/group/SOL/guides/qpopt.pdf)).

### Obtaining a initial feasible point

An initial feasible point can be obtained e.g. by performing `Phase-I` `Simplex` on the polyhedron `Ax ≤ b`:
```
using JuMP, Gurobi
# Choose Gurobi's primal simplex method
model = Model(solver=GurobiSolver(Presolve=0, Method=0))
@variable(model, x[1:size(A, 2)])
@constraint(model, A*x - b .<=0)
status = JuMP.solve(model)

x_init = getvalue(x)  # Initial point to be passed to our solver
```