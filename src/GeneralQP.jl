__precompile__(true)

module GeneralQP

using LinearAlgebra

include("linear_algebra.jl")
include("change_constraints.jl")
include("qp.jl")
export solve, UpdatableQR, NullspaceHessianLDL

end