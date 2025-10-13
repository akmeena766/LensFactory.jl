module ExternalEffects

# Inbuilt packages to use

# LensFactory modules to import
using ..Constants

# Functions to export
export potential! 
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: RV
   ψ_up = ψ + 0.5 * (kv + g1) * θx^2 + 0.5 * (kv - g1) * θy^2 + g2 * θx * θy
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: ROA
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = ψ[i, j] + 0.5 * (kv + g1) * θx[i, j]^2 + 0.5 * (kv - g1) * θy[i, j]^2 + g2 * θx[i, j] * θy[i, j]
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: RV
   ψx_up = ψx + (kv + g1) * θx + g2 * θy
   ψy_up = ψy + (kv - g1) * θy + g2 * θx
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: ROA
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
      ψx[i, j] = ψx[i, j] + (kv + g1) * θx[i, j] + g2 * θy[i, j]
      ψy[i, j] = ψy[i, j] + (kv - g1) * θy[i, j] + g2 * θx[i, j]
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: RV
   ψxx_up = ψxx + (kv + g1)
   ψyy_up = ψyy + (kv - g1)
   ψxy_up = ψxy + g2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kv::RV, g1::RV, g2::RV) where T <: ROA
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψxx[i, j] = ψxx[i, j] + (kv + g1)
         ψyy[i, j] = ψyy[i, j] + (kv - g1)
         ψxy[i, j] = ψxy[i, j] + g2
      end
   end
end

end