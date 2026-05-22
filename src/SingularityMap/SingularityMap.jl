module SingularityMap


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Lenses


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export from_lens
export from_jacobian


# --------------------------------------------------------------------------------------------------
# Local functions to include
# --------------------------------------------------------------------------------------------------
include("./eigen.jl")
include("./extrema.jl")
include("./a3lines.jl")


function from_jacobian(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T) where T <: Matrix{Float64}
   # Dimensions
   nx, ny = size(θx)

   # Construct a single 3D array for jacobian
   ddf = zeros(Float64, nx, ny, 3)
   ddf[:, :, 1] = ψxx
   ddf[:, :, 2] = ψyy
   ddf[:, :, 3] = ψxy
   # println(findall(isinf, ddf))
   # println(minimum(ddf), maximum(ddf))
   # # Get eigenvalues and eigenvectors
   e1 = zeros(Float64, nx, ny)
   e2 = zeros(Float64, nx, ny)
   q1 = zeros(Float64, nx, ny, 2)
   q2 = zeros(Float64, nx, ny, 2)
   _eigen!(e1, e2, q1, q2, ddf)

   # Get gradient of eigenvalues
   de1 = zeros(Float64, nx, ny, 2)
   de2 = zeros(Float64, nx, ny, 2)
   _gradeigen!(de1, de2, e1, e2)

   # Get A3-lines
   na31::Int64 = 0
   na32::Int64 = 0
   ra31 = zeros(Float64, round(Int64, nx*ny/8), 2)
   ra32 = zeros(Float64, round(Int64, nx*ny/8), 2)
   na31 = _a3lines!(na31, ra31, e1, de1, q1)
   na32 = _a3lines!(na32, ra32, e2, de2, q2)
   @show na31, na32
   return na31, na32, ra31, ra32


   # # Get extrema, minima, and maxima points
   # nex1::Int64 = 0
   # nex2::Int64 = 0
   # rex1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # rex2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # nex1 = _extrema!(nex1, rex1, de1)
   # nex2 = _extrema!(nex2, rex2, de2)
   # println("Extrema: ", nex1, " ", nex2)

   # nmin1::Int64 = 0
   # nmin2::Int64 = 0
   # rmin1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # rmin2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # _minima!(nmin1, rmin1, e1, nex1, rex1)
   # _minima!(nmin2, rmin2, e2, nex2, rex2)
   # println("Minima: ", nex1, " ", nex2)

   # nmax1::Int64 = 0
   # nmax2::Int64 = 0
   # rmax1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # rmax2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # _maxima!(nmax1, rmax1, e1, nex1, rex1)
   # _maxima!(nmax2, rmax2, e2, nex2, rex2)
   # println("Maxima: ", nex1, " ", nex2)
   # println(rmax1[1, :])
end

function from_lens(lens::Lenses.AbstractLens, θx::T, θy::T) where T <: Matrix{Float64}
   # Get jacobian using lens
   ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, θx, θy)

   # Pass jacobian map to from_jacobian function
   return from_jacobian(ψxx, ψyy, ψxy, θx, θy)
end


end