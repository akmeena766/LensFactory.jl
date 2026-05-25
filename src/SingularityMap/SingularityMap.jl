module SingularityMap


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Lenses
using ..LFUtils


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


# --------------------------------------------------------------------------------------------------
# This is the common function that will be used by both from_jacobian and from_lens
# --------------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   # Remove points with a_crit > adis and covert positions from pixel to arcsecond
   # -----------------------------------------------------------------------------------------------
   na31_ref::Int64 = 0
   ra31_ref = zeros(Float64, round(Int64, nx * ny / 8), 2)
   for i in 1:na31
      e1_val = PolygonOps.bilinear_interpolation(ra31[i, 1], ra31[i, 2], e1)
      if adis * e1_val ≥ 1.0
         na31_ref = na31_ref + 1

         # Convert from pixel to arcsecond
         θx_value = PolygonOps.bilinear_interpolation(ra31[i, 1], ra31[i, 2], θx)
         θy_value = PolygonOps.bilinear_interpolation(ra31[i, 1], ra31[i, 2], θy)
         ra31_ref[na31_ref, 1] = θx_value
         ra31_ref[na31_ref, 2] = θy_value
      end
   end

   na32_ref::Int64 = 0
   ra32_ref = zeros(Float64, round(Int64, nx * ny / 8), 2)
   for i in 1:na32
      e2_val = PolygonOps.bilinear_interpolation(ra32[i, 1], ra32[i, 2], e2)
      if adis * e2_val ≥ 1.0
         na32_ref = na32_ref + 1

         θx_value = PolygonOps.bilinear_interpolation(ra32[i, 1], ra32[i, 2], θx)
         θy_value = PolygonOps.bilinear_interpolation(ra32[i, 1], ra32[i, 2], θy)
         ra32_ref[na32_ref, 1] = θx_value
         ra32_ref[na32_ref, 2] = θy_value
      end
   end


   nd4_ref::Int64 = 0
   rd4_ref = zeros(Float64, round(Int64, nx * ny / 8), 2)
   for i in 1:nd4
      e1_val = PolygonOps.bilinear_interpolation(rd4[i, 1], rd4[i, 2], e1)
      e2_val = PolygonOps.bilinear_interpolation(rd4[i, 1], rd4[i, 2], e2)
      if adis * e1_val ≥ 1.0 && adis * e2_val ≥ 1.0
         nd4_ref = nd4_ref + 1

         θx_value = PolygonOps.bilinear_interpolation(rd4[i, 1], rd4[i, 2], θx)
         θy_value = PolygonOps.bilinear_interpolation(rd4[i, 1], rd4[i, 2], θy)
         rd4_ref[nd4_ref, 1] = θx_value
         rd4_ref[nd4_ref, 2] = θy_value
      end
   end
   
   return ra31_ref[1:na31_ref, :], ra32_ref[1:na32_ref, :], rd4_ref[1:nd4_ref, :]
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
   ra31, ra32, rd4 = _common(ddf, θx, θy; adis=adis, buffer=buffer)

   return ra31, ra32, rd4
end


end