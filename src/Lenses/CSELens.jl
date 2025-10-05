module CSELens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   
end

end