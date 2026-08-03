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
   return log(θx^2 + θy^2) / (2π)
end


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T, nl::Int64) where {U<:Real, S<:Real, T<:Vector{<:Real}}
   ψ_up = ψ
   @inbounds for k in 1:nl
      xp = θx - (θxc[k] + θpix[k]/2)
      xm = θx - (θxc[k] - θpix[k]/2)
      yp = θy - (θyc[k] + θpix[k]/2)
      ym = θy - (θyc[k] - θpix[k]/2)

      ψ_up = ψ_up + κ[k] *( _ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp) )
   end
   return ψ_up
end

function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T, nl::Int64) where {U<:ROA, S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for k in 1:nl
      @inbounds for j in ax2
         @inbounds for i in ax1
            xp = θx[i, j] - (θxc[k] + θpix[k]/2)
            xm = θx[i, j] - (θxc[k] - θpix[k]/2)
            yp = θy[i, j] - (θyc[k] + θpix[k]/2)
            ym = θy[i, j] - (θyc[k] - θpix[k]/2)
            
            ψ[i, j] = ψ[i, j] + κ[k] * (_ψ̃(xp, yp) + _ψ̃(xm, ym) - _ψ̃(xp, ym) - _ψ̃(xm, yp))
         end
      end
   end
   return nothing
end


function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T, nl::Int64) where {U<:Real, S<:Real, T<:Vector{<:Real}}
   ψx_up = ψx
   ψy_up = ψy
   @inbounds for k in 1:nl
      xp = θx - (θxc[k] + θpix[k]/2)
      xm = θx - (θxc[k] - θpix[k]/2)
      yp = θy - (θyc[k] + θpix[k]/2)
      ym = θy - (θyc[k] - θpix[k]/2)

      ψx_up = ψx_up + κ[k] * (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
      ψy_up = ψy_up + κ[k] * (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
   end
   return ψx_up, ψy_up
end

function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T, nl::Int64) where {U<:ROA, S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for k in 1:nl
      @inbounds for j in ax2
         @inbounds for i in ax1
            xp = θx[i, j] - (θxc[k] + θpix[k]/2)
            xm = θx[i, j] - (θxc[k] - θpix[k]/2)
            yp = θy[i, j] - (θyc[k] + θpix[k]/2)
            ym = θy[i, j] - (θyc[k] - θpix[k]/2)
            
            ψx[i, j] = ψx[i, j] + κ[k] * (_ψ̃x(xp, yp) + _ψ̃x(xm, ym) - _ψ̃x(xp, ym) - _ψ̃x(xm, yp))
            ψy[i, j] = ψy[i, j] + κ[k] * (_ψ̃y(xp, yp) + _ψ̃y(xm, ym) - _ψ̃y(xp, ym) - _ψ̃y(xm, yp))
         end
      end
   end
   return nothing
end


function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T, nl::Int64) where {U<:Real, S<:Real, T<:Vector{<:Real}}
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy
   @inbounds for k in 1:nl
      xp = θx - (θxc[k] + θpix[k]/2)
      xm = θx - (θxc[k] - θpix[k]/2)
      yp = θy - (θyc[k] + θpix[k]/2)
      ym = θy - (θyc[k] - θpix[k]/2)

      ψxx_up = ψxx_up + κ[k] * (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
      ψyy_up = ψyy_up + κ[k] * (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
      ψxy_up = ψxy_up + κ[k] * (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
   end
   return ψxx_up, ψyy_up, ψxy_up
end

function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, κ::T, θpix::T, nl::Int64) where {U<:ROA, S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for k in 1:nl
      @inbounds for j in ax2
         @inbounds for i in ax1
            xp = θx[i, j] - (θxc[k] + θpix[k]/2)
            xm = θx[i, j] - (θxc[k] - θpix[k]/2)
            yp = θy[i, j] - (θyc[k] + θpix[k]/2)
            ym = θy[i, j] - (θyc[k] - θpix[k]/2)
            
            ψxx[i, j] = ψxx[i, j] + κ[k] * (_ψ̃xx(xp, yp) + _ψ̃xx(xm, ym) - _ψ̃xx(xp, ym) - _ψ̃xx(xm, yp))
            ψyy[i, j] = ψyy[i, j] + κ[k] * (_ψ̃yy(xp, yp) + _ψ̃yy(xm, ym) - _ψ̃yy(xp, ym) - _ψ̃yy(xm, yp))
            ψxy[i, j] = ψxy[i, j] + κ[k] * (_ψ̃xy(xp, yp) + _ψ̃xy(xm, ym) - _ψ̃xy(xp, ym) - _ψ̃xy(xm, yp))
         end
      end
   end
   return nothing
end

end