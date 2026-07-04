module LieAlgebra

using LinearAlgebra
import LinearAlgebra: rank

include("type.jl")
export MatrixBasis

include("inner_products.jl")
export AbstractInnerProduct, inner_product, gram_matrix
export TraceInnerProduct

include("basis.jl")
export elements, matrix_size, basis_matrix, is_independent, assert_independent, coefficients, reconstruct, project, residual, in_span, is_closed, assert_closed

include("brackets.jl")
export commutator, bracket, bracket_coefficients

include("structure_constants.jl")
export structure_constants

include("adjoint.jl")
export adjoint_matrix, adjoint_matrices

include("killing.jl")
export killing_form

end
