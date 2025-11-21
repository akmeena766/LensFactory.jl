module ExternalEffects

# Inbuilt packages to use

# LensFactory modules to import
using ..Constants

# Functions to export
export potential! 
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: RV
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   ψ_up = ψ + 0.5 * (kappa + gamma1) * θx^2 + 0.5 * (kappa - gamma1) * θy^2 + gamma2 * θx * θy
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: ROA
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = ψ[i, j] + 0.5 * (kappa + gamma1) * θx[i, j]^2 + 0.5 * (kappa - gamma1) * θy[i, j]^2 + gamma2 * θx[i, j] * θy[i, j]
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: RV
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   ψx_up = ψx + (kappa + gamma1) * θx + gamma2 * θy
   ψy_up = ψy + (kappa - gamma1) * θy + gamma2 * θx
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: ROA
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
      ψx[i, j] = ψx[i, j] + (kappa + gamma1) * θx[i, j] + gamma2 * θy[i, j]
      ψy[i, j] = ψy[i, j] + (kappa - gamma1) * θy[i, j] + gamma2 * θx[i, j]
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: RV
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   ψxx_up = ψxx + (kappa + gamma1)
   ψyy_up = ψyy + (kappa - gamma1)
   ψxy_up = ψxy + gamma2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, kappa::RV, gamma::RV, angle::RV) where T <: ROA
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψxx[i, j] = ψxx[i, j] + (kappa + gamma1)
         ψyy[i, j] = ψyy[i, j] + (kappa - gamma1)
         ψxy[i, j] = ψxy[i, j] + gamma2
      end
   end
end

end