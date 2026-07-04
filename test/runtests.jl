using LieAlgebra
using Test

using LinearAlgebra
using Random

px = [0 1; 1 0]
py = [0 -im; im 0]
pz = [1 0; 0 -1]
I2 = Matrix(I, 2, 2)

Sx = kron(px, I2)
Sy = kron(py, I2)
Sz = kron(pz, I2)
Tx = kron(I2, px)
Ty = kron(I2, py)
Tz = kron(I2, pz)

include("matrix_basis.jl")
include("structure_constants.jl")
include("commutant.jl")
include("decomposition.jl")
include("orthogonalization.jl")
