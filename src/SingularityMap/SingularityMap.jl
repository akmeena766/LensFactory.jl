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
include("./a4points.jl")
include("./d4points.jl")


function _common(ddf::Array{Float64, 3}, θx::T, θy::T; adis::Float64=1.0, buffer::Int64=10) where T<:Matrix{Float64}
   # Dimensions
   nx, ny = size(θx)

   # Get eigenvalues and eigenvectors
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
   ra31 = zeros(Float64, round(Int64, nx * ny / 8), 2)
   ra32 = zeros(Float64, round(Int64, nx * ny / 8), 2)
   na31 = _a3lines!(na31, ra31, e1, de1, q1)
   na32 = _a3lines!(na32, ra32, e2, de2, q2)

   # Get extrema
   nex1::Int64 = 0
   nex2::Int64 = 0
   # rex1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # rex2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # nex1 = _extrema!(nex1, rex1, de1)
   # nex2 = _extrema!(nex2, rex2, de2)

   # Get minima
   nmin1::Int64 = 0
   nmin2::Int64 = 0
   # rmin1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # rmin2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # _minima!(nmin1, rmin1, e1, nex1, rex1)
   # _minima!(nmin2, rmin2, e2, nex2, rex2)

   # Get maxima
   nmax1::Int64 = 0
   nmax2::Int64 = 0
   # rmax1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # rmax2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   # _maxima!(nmax1, rmax1, e1, nex1, rex1)
   # _maxima!(nmax2, rmax2, e2, nex2, rex2)

   # Get D4-points
   nd4::Int64 = 0
   rd4 = zeros(Float64, round(Int64, nx * ny / 8), 2)
   nd4 = _d4points(nd4, rd4, ddf)
   
   return ra31[1:na31, :], ra32[1:na32, :], rd4[1:nd4, :]
end


"""
    from_jacobian(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T;
                  adis::Float64 = 1.0,
                  buffer::Int64 = 10) where T<:Matrix{Float64}
"""
function from_jacobian(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T; adis::Float64=1.0, buffer::Int64=10) where T<:Matrix{Float64}
   # Dimensions
   nx, ny = size(θx)

   # Construct a single 3D array for jacobian
   ddf = zeros(Float64, nx, ny, 3)
   ddf[:, :, 1] = ψxx
   ddf[:, :, 2] = ψyy
   ddf[:, :, 3] = ψxy

   # Get singularity map
   ra31, ra32, rd4 = _common(ddf, θx, θy; adis=adis, buffer=buffer)

   return ra31, ra32, rd4
end


"""
    from_lens(lens::Lenses.AbstractLens, θx::T, θy::T; 
              adis:: Float64      = 1.0,
              buffer::Int64       = 10,
              adaptive::Bool      = false,
              resolution::Float64 = 0.001) where T<:Matrix{Float64})
"""
function from_lens(lens::Lenses.AbstractLens, θx::T, θy::T;
                   adis::Float64       = 1.0,
                   buffer::Int64       = 10,
                   adaptive::Bool      = false,
                   resolution::Float64 = 0.001) where T<:Matrix{Float64}
   # Dimensions
   nx, ny = size(θx)
   pixel_size = θx[1, 2] - θx[1, 1]

   # Get jacobian using lens
   ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, θx, θy)
   
   # Combined array
   ddf = zeros(Float64, nx, ny, 3)
   ddf[:, :, 1] = ψxx
   ddf[:, :, 2] = ψyy
   ddf[:, :, 3] = ψxy

   


   # Pass jacobian map to from_jacobian function
   return _common(ψxx, ψyy, ψxy, θx, θy; )
end


end