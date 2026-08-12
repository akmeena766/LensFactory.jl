module PixelLens


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
   return log(θx^2 + θy^2) / (2π)
end


function _ψ̃_kernel(θx::T, θy::T, θxc::T, θyc::T, θpix::T) where T <: Real
   xp = θx - (θxc + θpix/2)
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)
   ym = θy - (θyc - θpix/2)

   return (_ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp))
end

function _ψ̃x_kernel(θx::T, θy::T, θxc::T, θyc::T, θpix::T) where T <: Real
   xp = θx - (θxc + θpix/2)
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)
   ym = θy - (θyc - θpix/2)

   return (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
end

function _ψ̃y_kernel(θx::T, θy::T, θxc::T, θyc::T, θpix::T) where T <: Real
   xp = θx - (θxc + θpix/2)
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)
   ym = θy - (θyc - θpix/2)

   return (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
end

function _ψ̃xx_kernel(θx::T, θy::T, θxc::T, θyc::T, θpix::T) where T <: Real
   xp = θx - (θxc + θpix/2)
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)
   ym = θy - (θyc - θpix/2)

   return (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
end

function _ψ̃yy_kernel(θx::T, θy::T, θxc::T, θyc::T, θpix::T) where T <: Real
   xp = θx - (θxc + θpix/2)
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)
   ym = θy - (θyc - θpix/2)

   return (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
end

function _ψ̃xy_kernel(θx::T, θy::T, θxc::T, θyc::T, θpix::T) where T <: Real
   xp = θx - (θxc + θpix/2)
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)
   ym = θy - (θyc - θpix/2)

   return (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
end


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
"""
    potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:Real, S<:Real, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:Real, S<:Real, T<:Real}
   xp = θx - (θxc + θpix/2)  
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)  
   ym = θy - (θyc - θpix/2)

   ψ_up = ψ + κ *( _ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp) )
   return ψ_up
end

"""
    potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate potential at given coordinates for a pixel lens and update the potential (ψ) in place.

# Arguments
- `ψ`   : Potential at given coordinates
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc` : x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc` : y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `κ`   : Convergence value (assuming ``D_{ds}/D_s = 1``).
- `θpix`: Pixel size (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the potential (ψ) in place.

"""
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:ROA, S<:ROA, T<:Real}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         xp = θx[i, j] - (θxc + θpix/2)
         xm = θx[i, j] - (θxc - θpix/2)
         yp = θy[i, j] - (θyc + θpix/2)
         ym = θy[i, j] - (θyc - θpix/2)
         
         ψ[i, j] = ψ[i, j] + κ * (_ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp))
      end
   end
   return nothing
end


"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:Real, S<:Real, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:Real, S<:Real, T<:Real}
   xp = θx - (θxc + θpix/2)  
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)  
   ym = θy - (θyc - θpix/2)

   ψx_up = ψx + κ * (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
   ψy_up = ψy + κ * (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate deflection at given coordinates for a pixel lens and update the deflection components
(ψx, ψy) in place.

# Arguments
- `ψx`  : x-component of the deflection at given coordinates
- `ψy`  : y-component of the deflection at given coordinates
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc` : x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc` : y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `κ`   : Convergence value (assuming ``D_{ds}/D_s = 1``).
- `θpix`: Pixel size (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the deflection (ψx, ψy) in place.

"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:ROA, S<:ROA, T<:Real}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         xp = θx[i, j] - (θxc + θpix/2)
         xm = θx[i, j] - (θxc - θpix/2)
         yp = θy[i, j] - (θyc + θpix/2)
         ym = θy[i, j] - (θyc - θpix/2)
         
         ψx[i, j] = ψx[i, j] + κ * (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
         ψy[i, j] = ψy[i, j] + κ * (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:Real, S<:Real, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:Real, S<:Real, T<:Real}
   xp = θx - (θxc + θpix/2)  
   xm = θx - (θxc - θpix/2)
   yp = θy - (θyc + θpix/2)  
   ym = θy - (θyc - θpix/2)

   ψxx_up = ψxx + κ * (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
   ψyy_up = ψyy + κ * (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
   ψxy_up = ψxy + κ * (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate jacobian at given coordinates for a pixel lens and update the jacobian components 
(ψxx, ψyy, ψxy) in place.

# Arguments
- `ψxx` : x-component of the jacobian at given coordinates
- `ψyy` : y-component of the jacobian at given coordinates
- `ψxy` : xy-component of the jacobian at given coordinates
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc` : x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc` : y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `κ`   : Convergence value (assuming ``D_{ds}/D_s = 1``).
- `θpix`: Pixel size (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the jacobian (ψxx, ψyy, ψxy) in place.
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T) where {U<:ROA, S<:ROA, T<:Real}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         xp = θx[i, j] - (θxc + θpix/2)
         xm = θx[i, j] - (θxc - θpix/2)
         yp = θy[i, j] - (θyc + θpix/2)
         ym = θy[i, j] - (θyc - θpix/2)
         
         ψxx[i, j] = ψxx[i, j] + κ * (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
         ψyy[i, j] = ψyy[i, j] + κ * (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
         ψxy[i, j] = ψxy[i, j] + κ * (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
      end
   end
   return nothing
end

end