using LieAlgebra
using Test

using LinearAlgebra
using Random

px = [0 1; 1 0]
py = [0 -im; im 0]
pz = [1 0; 0 -1]
I2 = Matrix(I, 2, 2)

include("matrix_basis.jl")
include("structure_constants.jl")
include("commutant.jl")
include("decomposition.jl")
