__precompile__(true)

module GeneralQP

using LinearAlgebra
using Polynomials

include("linear_algebra.jl")
include("change_constraints.jl")
include("printing.jl")
include("qp.jl")
export solve
export UpdatableQR, NullspaceHessianLDL
export remove_constraint!, add_constraint!

end