module MultiPixelLens


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export potential! 
export deflection!
export jacobian!


# --------------------------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------------------------
@inline function _atan(θx::T, θy::T) where T <: Real
   if abs(θx) < eps(T)
      return sign(θy) * π/2
   end
   return atan(θy / θx)
end

@inline function _ψ̃(θx::T, θy::T) where T <: Real
   θ2  = θx^2 + θy^2
   if θ2 < eps(T)
      return zero(T)
   end
   return (θx^2 *_atan(θx, θy) + θy^2 * _atan(θy, θx) + θx * θy * log(θ2) - 3.0 * θx * θy) / (2π)
end

@inline function _ψ̃x(θx::T, θy::T) where T <: Real
   θ2  = θx^2 + θy^2
   if θ2 < eps(T)
      return (θx * _atan(θx, θy) - θy) / π
   end
   return (θx * _atan(θx, θy) + θy * log(θ2) / 2 - θy) / π
end

@inline function _ψ̃y(θx::T, θy::T) where T <: Real
   θ2 = θx^2 + θy^2
   if θ2 < eps(T)
      return (θy * _atan(θy, θx) - θx) / π
   end
   return (θy * _atan(θy, θx) + θx * log(θ2) / 2 - θx) / π
end

@inline function _ψ̃xx(θx::T, θy::T) where T <: Real
   return _atan(θx, θy) / π
end

@inline function _ψ̃yy(θx::T, θy::T) where T <: Real
   return _atan(θy, θx) / π
end

@inline function _ψ̃xy(θx::T, θy::T) where T <: Real
   θ2 = θx^2 + θy^2
   if θ2 < eps(T)
      return zero(T)
   end
   return log(θ2) / (2π)
end

function _check_grid(θxc::AbstractMatrix{<:Real}, θyc::AbstractMatrix{<:Real}, θpix::Real)
   nx, ny = size(θxc)
   if size(θyc) != (nx, ny)
      throw(ArgumentError("x_c and y_c must have the same size; got $(size(θxc)) and $(size(θyc))."))
   end
   if nx < 1 || ny < 1 || !(θpix > 0)
      throw(ArgumentError("The grid needs at least one pixel and a positive pixel size."))
   end
   tol = sqrt(eps(Float64)) * θpix
   for j in 1:ny, i in 1:nx
      if abs(θxc[i, j] - θxc[i, 1]) > tol || abs(θyc[i, j] - θyc[1, j]) > tol
         throw(ArgumentError("x_c must vary only along the first dimension and y_c only along the " * 
                             "second (as from get_meshgrid)."))
      end
   end
   for i in 2:nx
      if abs(θxc[i, 1] - θxc[i - 1, 1] - θpix) > tol
         throw(ArgumentError("Pixel centres must be spaced by the pixel size along x."))
      end
   end
   for j in 2:ny
      if abs(θyc[1, j] - θyc[1, j - 1] - θpix) > tol
         throw(ArgumentError("Pixel centres must be spaced by the pixel size along y."))
      end
   end
   return nothing
end

function _corners(θxc::AbstractMatrix{<:Real}, θyc::AbstractMatrix{<:Real}, θpix::Real)
   nx, ny = size(θxc)
   h  = θpix / 2
   cx = Vector{Float64}(undef, nx + 1)
   cy = Vector{Float64}(undef, ny + 1)
   @inbounds for a in 1:nx
      cx[a] = θxc[a, 1] - h            # left edge of pixel column a
   end
   cx[nx + 1] = θxc[nx, 1] + h         # right edge of the last column
 
   @inbounds for b in 1:ny
      cy[b] = θyc[1, b] - h            # bottom edge of pixel row b
   end
   cy[ny + 1] = θyc[1, ny] + h         # top edge of the last row
   return cx, cy
end


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::Real) where {U<:Real, S<:Real, T<:Matrix{<:Real}}
   # Get size of the multi-pixel lens
   nx, ny = size(κ)

   # Get corners
   cx, cy = _corners(θxc, θyc, θpix)

   ψ_up = ψ
   @inbounds for j1 in 1:ny
      @inbounds for i1 in 1:nx
         xp = θx - cx[i1 + 1]
         xm = θx - cx[i1]
         yp = θy - cy[j1 + 1]
         ym = θy - cy[j1]

         ψ_up = ψ_up + κ[i1, j1] *( _ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp) )
      end
   end
   return ψ_up
end

function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::Real) where {U<:ROA, S<:ROA, T<:Matrix{<:Real}}
   # Get size of the multi-pixel lens
   nx, ny = size(κ)

   # Get corners
   cx, cy = _corners(θxc, θyc, θpix)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   Threads.@threads for j in ax2
      @inbounds for i in ax1
         @inbounds for j1 in 1:ny
            @inbounds for i1 in 1:nx
               xp = θx[i, j] - cx[i1 + 1]
               xm = θx[i, j] - cx[i1]
               yp = θy[i, j] - cy[j1 + 1]
               ym = θy[i, j] - cy[j1]

               ψ[i, j] = ψ[i, j] + κ[i1, j1] * (_ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp))
            end
         end
      end
   end
   return nothing
end


function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::Real) where {U<:Real, S<:Real, T<:Matrix{<:Real}}
   # Get size of the multi-pixel lens
   nx, ny = size(κ)

   # Get corners
   cx, cy = _corners(θxc, θyc, θpix)
   
   ψx_up = ψx
   ψy_up = ψy
   @inbounds for j1 in 1:ny
      @inbounds for i1 in 1:nx
         xp = θx - cx[i1 + 1]
         xm = θx - cx[i1]
         yp = θy - cy[j1 + 1]
         ym = θy - cy[j1]

         ψx_up = ψx_up + κ[i1, j1] * (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
         ψy_up = ψy_up + κ[i1, j1] * (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
      end
   end
   return ψx_up, ψy_up
end

function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::Real) where {U<:ROA, S<:ROA, T<:Matrix{<:Real}}
   # Get size of the multi-pixel lens
   nx, ny = size(κ)

   # Get corners
   cx, cy = _corners(θxc, θyc, θpix)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   Threads.@threads for j in ax2
      @inbounds for i in ax1
         @inbounds for j1 in 1:ny
            @inbounds for i1 in 1:nx
               xp = θx[i, j] - cx[i1 + 1]
               xm = θx[i, j] - cx[i1]
               yp = θy[i, j] - cy[j1 + 1]
               ym = θy[i, j] - cy[j1]

               ψx[i, j] = ψx[i, j] + κ[i1, j1] * (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
               ψy[i, j] = ψy[i, j] + κ[i1, j1] * (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
            end
         end
      end
   end
   return nothing
end


function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::Real) where {U<:Real, S<:Real, T<:Matrix{<:Real}}
   # Get size of the multi-pixel lens
   nx, ny = size(κ)

   # Get corners
   cx, cy = _corners(θxc, θyc, θpix)
   
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy
   @inbounds for j1 in 1:ny
      @inbounds for i1 in 1:nx
         xp = θx - cx[i1 + 1]
         xm = θx - cx[i1]
         yp = θy - cy[j1 + 1]
         ym = θy - cy[j1]

         ψxx_up = ψxx_up + κ[i1, j1] * (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
         ψyy_up = ψyy_up + κ[i1, j1] * (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
         ψxy_up = ψxy_up + κ[i1, j1] * (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
      end
   end
   return ψxx_up, ψyy_up, ψxy_up
end

function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::Real) where {U<:ROA, S<:ROA, T<:Matrix{<:Real}}
   # Get size of the multi-pixel lens
   nx, ny = size(κ)

   # Get corners
   cx, cy = _corners(θxc, θyc, θpix)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   Threads.@threads for j in ax2
      @inbounds for i in ax1
         @inbounds for j1 in 1:ny
            @inbounds for i1 in 1:nx
               xp = θx[i, j] - cx[i1 + 1]
               xm = θx[i, j] - cx[i1]
               yp = θy[i, j] - cy[j1 + 1]
               ym = θy[i, j] - cy[j1]
            
               ψxx[i, j] = ψxx[i, j] + κ[i1, j1] * (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
               ψyy[i, j] = ψyy[i, j] + κ[i1, j1] * (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
               ψxy[i, j] = ψxy[i, j] + κ[i1, j1] * (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
            end
         end
      end
   end
   return nothing
end


# --------------------------------------------------------------------------------------------------
# Internal functions relevant for lens model fitting
# --------------------------------------------------------------------------------------------------
# Fill column p of K with the four-corner differences of the corner function f, evaluated once on
# the (nx+1) × (ny+1) lattice of pixel corners. Neighbouring pixels share corners, so this needs
# about 4× fewer evaluations of f than evaluating four corners per pixel. Two lattice rows are kept.
@inline _kernel_column!(::Nothing, f::F, p, u, v, prev, curr, nx, ny) where F = nothing
 
function _kernel_column!(K::Matrix{<:Real}, f::F, p::Int64, u::Vector{Float64}, v::Vector{Float64}, 
                         prev::Vector{Float64}, curr::Vector{Float64}, nx::Int64, ny::Int64) where F
   @inbounds for a in 1:(nx + 1)
      prev[a] = f(u[a], v[1])
   end
   @inbounds for j in 1:ny
      vj = v[j + 1]
      for a in 1:(nx + 1)
         curr[a] = f(u[a], vj)
      end
      off = (j - 1) * nx
      for i in 1:nx
         K[off + i, p] = curr[i + 1] + prev[i] - prev[i + 1] - curr[i]
      end
      prev, curr = curr, prev
   end
   return nothing
end


function _kernel!(Kψ::Union{Nothing, Matrix{<:Real}}, 
                  Kαx::Union{Nothing, Matrix{<:Real}}, 
                  Kαy::Union{Nothing, Matrix{<:Real}}, 
                  Kψxx::Union{Nothing, Matrix{<:Real}}, 
                  Kψyy::Union{Nothing, Matrix{<:Real}}, 
                  Kψxy::Union{Nothing, Matrix{<:Real}},
                  θx::AbstractVector{<:Real}, 
                  θy::AbstractVector{<:Real}, 
                  θxc::AbstractVector{<:Real}, 
                  θyc::AbstractVector{<:Real}, 
                  θpix::Real)
   # Dimensions
   n      = length(θx)
   nx, ny = size(θxc)

   # Input checks
   if length(θy) != n
      throw(ArgumentError("θx and θy must have the same length."))
   end
   _check_grid(θxc, θyc, θpix)

   for K in (Kψ, Kαx, Kαy, Kψxx, Kψyy, Kψxy)
      if K !== nothing && size(K) != (nx * ny, n)
         throw(DimensionMismatch("Kernel matrices must be of size ($(nx * ny), $n); got $(size(K))."))
      end
   end

   # Pixel-edge coordinates along each axis (nx + 1 and ny + 1 of them)
   cx, cy = _corners(θxc, θyc, θpix)

   # Local arrays
   u    = Vector{Float64}(undef, nx + 1)
   v    = Vector{Float64}(undef, ny + 1)
   prev = Vector{Float64}(undef, nx + 1)
   curr = Vector{Float64}(undef, nx + 1)

   for p in 1:n
      @inbounds for a in 1:(nx + 1)
         u[a] = θx[p] - cx[a]
      end
      @inbounds for b in 1:(ny + 1)
         v[b] = θy[p] - cy[b]
      end
      _kernel_column!(Kψ,   _ψ̃,   p, u, v, prev, curr, nx, ny)
      _kernel_column!(Kαx,  _ψ̃x,  p, u, v, prev, curr, nx, ny)
      _kernel_column!(Kαy,  _ψ̃y,  p, u, v, prev, curr, nx, ny)
      _kernel_column!(Kψxx, _ψ̃xx, p, u, v, prev, curr, nx, ny)
      _kernel_column!(Kψyy, _ψ̃yy, p, u, v, prev, curr, nx, ny)
      _kernel_column!(Kψxy, _ψ̃xy, p, u, v, prev, curr, nx, ny)
   end
   return nothing
end

end