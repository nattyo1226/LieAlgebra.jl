module LieAlgebra

using LinearAlgebra
import LinearAlgebra: rank
using Random

include("type.jl")
export MatrixBasis, Decomposition

include("inner_products.jl")
export AbstractInnerProduct, inner_product, gram_matrix
export TraceInnerProduct

include("basis.jl")
export elements, matrix_size, basis_matrix, is_independent, assert_independent, coefficients, reconstruct, project, residual, in_span, is_closed, assert_closed

include("brackets.jl")
export commutator, lie_bracket, bracket, bracket_coefficients

include("structure_constants.jl")
export structure_constants

include("adjoint.jl")
export adjoint_matrix, adjoint_matrices

include("killing.jl")
export killing_form

include("commutant.jl")
export commutant

include("decomposition.jl")
export generic_commutant_element, decompose, transform_basis, transformed_basis, component_basis, component_bases

end
