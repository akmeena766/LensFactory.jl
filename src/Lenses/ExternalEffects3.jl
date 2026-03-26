module ExternalEffects3


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


"""
    potential!(ψ::T, θx::T, θy::T, δ::RV, θ::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, δ::RV, θ::RV) where T <: RV

end


"""
    potential!(ψ::T, θx::T, θy::T, δ::RV, θ::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, δ::RV, θ::RV) where T <: ROA
   
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: RV
   
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: ROA
   
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: RV
   
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, θ::RV) where T <: ROA
   
end

end