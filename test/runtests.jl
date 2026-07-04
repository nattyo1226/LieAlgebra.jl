using LieAlgebra
using Test

using LinearAlgebra

px = [0 1; 1 0]
py = [0 -im; im 0]
pz = [1 0; 0 -1]
B = MatrixBasis([px, py, pz])

include("matrix_basis.jl")
include("structure_constants.jl")
include("commutant.jl")
