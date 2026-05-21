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
export from_potential
export from_deflection
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

   # Get extrema, minima, and maxima points
   nex1::Int64 = 0
   nex2::Int64 = 0
   rex1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   rex2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   nex1 = _extrema!(nex1, rex1, de1)
   nex2 = _extrema!(nex2, rex2, de2)
   println("Extrema: ", nex1, " ", nex2)

   nmin1::Int64 = 0
   nmin2::Int64 = 0
   rmin1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   rmin2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   _minima!(nmin1, rmin1, e1, nex1, rex1)
   _minima!(nmin2, rmin2, e2, nex2, rex2)
   println("Minima: ", nex1, " ", nex2)

   nmax1::Int64 = 0
   nmax2::Int64 = 0
   rmax1 = zeros(Float64, round(Int64, nx*ny/8), 2)
   rmax2 = zeros(Float64, round(Int64, nx*ny/8), 2)
   _maxima!(nmax1, rmax1, e1, nex1, rex1)
   _maxima!(nmax2, rmax2, e2, nex2, rex2)
   println("Maxima: ", nex1, " ", nex2)
   println(rmax1[1, :])

   # Get A3-lines
   na31::Int64 = 0
   na32::Int64 = 0
   ra31 = zeros(Float64, round(Int64, nx*ny/8), 2)
   ra32 = zeros(Float64, round(Int64, nx*ny/8), 2)
   _a3lines!(na31, ra31, e1, de1, q1)
   _a3lines!(na32, ra32, e2, de2, q2)
   println(na31, na32)
end


function from_deflection(ψx::T, ψy::T, θx::T, θy::T) where T <: Matrix{Float64}
   # Pixel size along x- and y-axis
   dx = abs(θx[2, 1] - θx[1, 1])
   dy = abs(θy[1, 2] - θy[1, 1])

   # Calculate jacobian map using finite difference
   ψxx = zero(θx)
   ψyy = zero(θx)
   ψxy = zero(θx)

   ax1, ax2 = axes(θx, 1)[3:end-2], axes(θx, 2)[3:end-2]
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψxx[i, j] = (ψx[i+1, j] - ψx[i-1, j]) / (2.0 * dx)
         ψyy[i, j] = (ψy[i, j+1] - ψy[i, j-1]) / (2.0 * dy)
         ψxy[i, j] = 0.5 * ( (ψx[i, j+1] - ψx[i, j-1]) / (2.0 * dy) + (ψy[i+1, j] - ψy[i-1, j]) / (2.0 * dx) )
      end
   end

   # Pass jacobian map to from_jacobian function
   return from_jacobian(ψxx, ψyy, ψxy, θx, θy)
end


function from_potential(ψ::T, θx::T, θy::T) where T <: Matrix{Float64}
   # Pixel size along x- and y-axis
   dx = abs(θx[2, 1] - θx[1, 1])
   dy = abs(θy[1, 2] - θy[1, 1])

   # Calculate deflection map using finite difference
   ψx = zero(θx)
   ψy = zero(θx)

   ax1, ax2 = axes(θx, 1)[2:end-1], axes(θx, 2)[2:end-1]
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψx[i, j] = (ψ[i+1, j] - ψ[i-1, j]) / (2.0 * dx)
         ψy[i, j] = (ψ[i, j+1] - ψ[i, j-1]) / (2.0 * dy)
      end
   end

   # Pass deflection map to from_deflecion function and return 
   return return from_deflection(ψx, ψy, θx, θy)
end


function from_lens(lens::Lenses.AbstractLens, θx::T, θy::T) where T <: Matrix{Float64}
   # Get jacobian using lens
   ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, θx, θy)

   # Pass jacobian map to from_jacobian function
   return from_jacobian(ψxx, ψyy, ψxy, θx, θy)
end


end