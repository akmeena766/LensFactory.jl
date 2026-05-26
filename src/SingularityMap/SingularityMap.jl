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
   # -----------------------------------------------------------------------------------------------
   # Common qunatities
   # -----------------------------------------------------------------------------------------------
   # Dimensions
   nx, ny = size(θx)

   # Get eigenvalues and eigenvectors
   e1 = zeros(Float64, nx, ny)
   e2 = zeros(Float64, nx, ny)
   q1 = zeros(Float64, 2, nx, ny)
   q2 = zeros(Float64, 2, nx, ny)
   _eigen!(e1, e2, q1, q2, ddf)

   # Get gradient of eigenvalues
   de1 = zeros(Float64, 2, nx, ny)
   de2 = zeros(Float64, 2, nx, ny)
   _gradeigen!(de1, de2, e1, e2)

   # -----------------------------------------------------------------------------------------------
   # Get A3-lines
   # -----------------------------------------------------------------------------------------------
   na31::Int64 = 0
   na32::Int64 = 0
   ra31 = zeros(Float64, 2, round(Int64, nx*ny/8))
   ra32 = zeros(Float64, 2, round(Int64, nx*ny/8))
   na31 = _a3lines!(na31, ra31, e1, de1, q1)
   na32 = _a3lines!(na32, ra32, e2, de2, q2)

   # -----------------------------------------------------------------------------------------------
   # Get A4-points
   # -----------------------------------------------------------------------------------------------
   # Get extrema
   nex1::Int64 = 0
   nex2::Int64 = 0
   rex1 = zeros(Float64, 2, round(Int64, nx*ny/8))
   rex2 = zeros(Float64, 2, round(Int64, nx*ny/8))
   nex1 = _extrema!(nex1, rex1, de1)
   nex2 = _extrema!(nex2, rex2, de2)

   # Get minima
   nmin1::Int64 = 0
   nmin2::Int64 = 0
   rmin1 = zeros(Float64, 2, round(Int64, nx*ny/8))
   rmin2 = zeros(Float64, 2, round(Int64, nx*ny/8))
   nmin1 = _minima!(nmin1, rmin1, e1, nex1, rex1)
   nmin2 = _minima!(nmin2, rmin2, e2, nex2, rex2)

   # Get maxima
   nmax1::Int64 = 0
   nmax2::Int64 = 0
   rmax1 = zeros(Float64, 2, round(Int64, nx*ny/8))
   rmax2 = zeros(Float64, 2, round(Int64, nx*ny/8))
   nmax1 = _maxima!(nmax1, rmax1, e1, nex1, rex1)
   nmax2 = _maxima!(nmax2, rmax2, e2, nex2, rex2)

   # Get A4-points
   na41::Int64 = 0
   na42::Int64 = 0
   ra41 = zeros(Float64, 2, round(Int64, nx*ny/8))
   ra42 = zeros(Float64, 2, round(Int64, nx*ny/8))
   na41 = _a4points!(na41, ra41, e1, na31, ra31, nex1, rex1, nmax1, rmax1, nmin1, rmin1)
   # na42 = _a4points!(na42, ra42, e2, na32, ra32, nex2, rex2, nmax2, rmax2, nmin2, rmin2)

   # return ra31_ref[1:na31 :], ra32_ref[1:na32, :], ra41[1:na41, :], ra42[1:na42, :], rd4_ref[1:nd4, :]   

   # -----------------------------------------------------------------------------------------------
   # Get D4-points
   # -----------------------------------------------------------------------------------------------
   nd4::Int64 = 0
   rd4 = zeros(Float64, 2, round(Int64, nx*ny/8))
   nd4 = _d4points(nd4, rd4, ddf)

   # -----------------------------------------------------------------------------------------------
   # Remove points with a_crit > adis and covert positions from pixel to arcsecond
   # -----------------------------------------------------------------------------------------------
   na31_ref::Int64 = 0
   ra31_ref = zeros(Float64, 2, round(Int64, nx*ny/8))
   for i in 1:na31
      e1_val = PolygonOps.bilinear_interpolation(ra31[1, i], ra31[2, i], e1)
      if adis * e1_val ≥ 1.0
         na31_ref = na31_ref + 1

         # Convert from pixel to arcsecond
         θx_value = PolygonOps.bilinear_interpolation(ra31[1, i], ra31[2, i], θx)
         θy_value = PolygonOps.bilinear_interpolation(ra31[1, i], ra31[2, i], θy)
         ra31_ref[1, na31_ref] = θx_value
         ra31_ref[2, na31_ref] = θy_value
      end
   end

   na32_ref::Int64 = 0
   ra32_ref = zeros(Float64, 2, round(Int64, nx*ny/8))
   for i in 1:na32
      e2_val = PolygonOps.bilinear_interpolation(ra32[1, i], ra32[2, i], e2)
      if adis * e2_val ≥ 1.0
         na32_ref = na32_ref + 1

         θx_value = PolygonOps.bilinear_interpolation(ra32[1, i], ra32[2, i], θx)
         θy_value = PolygonOps.bilinear_interpolation(ra32[1, i], ra32[2, i], θy)
         ra32_ref[1, na32_ref] = θx_value
         ra32_ref[2, na32_ref] = θy_value
      end
   end


   nd4_ref::Int64 = 0
   rd4_ref = zeros(Float64, 2, round(Int64, nx*ny/8))
   for i in 1:nd4
      e1_val = PolygonOps.bilinear_interpolation(rd4[1, i], rd4[2, i], e1)
      e2_val = PolygonOps.bilinear_interpolation(rd4[1, i], rd4[2, i], e2)
      if adis * e1_val ≥ 1.0 && adis * e2_val ≥ 1.0
         nd4_ref = nd4_ref + 1

         θx_value = PolygonOps.bilinear_interpolation(rd4[1, i], rd4[2, i], θx)
         θy_value = PolygonOps.bilinear_interpolation(rd4[1, i], rd4[2, i], θy)
         rd4_ref[1, nd4_ref] = θx_value
         rd4_ref[2, nd4_ref] = θy_value
      end
   end
   
   return ra31_ref[:, 1:na31_ref], ra32_ref[:, 1:na32_ref], rd4_ref[:, 1:nd4_ref]
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
   ddf = zeros(Float64, 3, nx, ny)
   ddf[1, :, :] = ψxx
   ddf[2, :, :] = ψyy
   ddf[3, :, :] = ψxy

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
   ddf = zeros(Float64, 3, nx, ny)
   ddf[1, :, :] = ψxx
   ddf[2, :, :] = ψyy
   ddf[3, :, :] = ψxy

   # Pass jacobian map to from_jacobian function
   ra31, ra32, rd4 = _common(ddf, θx, θy; adis=adis, buffer=buffer)

   return ra31, ra32, rd4
end


end