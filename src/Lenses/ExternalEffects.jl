module ExternalEffects

# Inbuilt packages to use

# LensFactory modules to import
using ..Constants

# Functions to export
export potential! 
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: RV
   f1 = 0.5 * (kappa + gamma1)
   f2 = 0.5 * (kappa - gamma1)

   ψ_up = ψ + f1 * θx^2 + f2 * θy^2 + gamma2 * θx * θy
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: ROA
   f1 = 0.5 * (kappa + gamma1)
   f2 = 0.5 * (kappa - gamma1)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = ψ[i, j] + f1 * θx[i, j]^2 + f2 * θy[i, j]^2 + gamma2 * θx[i, j] * θy[i, j]
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: RV
   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ψx_up = ψx + f1 * θx + gamma2 * θy
   ψy_up = ψy + f2 * θy + gamma2 * θx
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: ROA
   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
      ψx[i, j] = ψx[i, j] + f1 * θx[i, j] + gamma2 * θy[i, j]
      ψy[i, j] = ψy[i, j] + f2 * θy[i, j] + gamma2 * θx[i, j]
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: RV
   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ψxx_up = ψxx + f1
   ψyy_up = ψyy + f2
   ψxy_up = ψxy + gamma2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma1::RV, gamma2::RV) where T <: ROA
   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψxx[i, j] = ψxx[i, j] + f1
         ψyy[i, j] = ψyy[i, j] + f2
         ψxy[i, j] = ψxy[i, j] + gamma2
      end
   end
end

end