"""
    Sources
"""
module Sources

# Julia packages to import
using SpecialFunctions

# LensFactory modules to use
using ..Constants

# Source profiles to export
export disk
export gaussian
export sersic

"""
    disk(θx::Matrix{<:RV}, θy::Matrix{<:RV}, θr::RV, β::NTuple{2, RV}; A::RV=1.0) --> Matrix{<:RV}

Creates a disk source profile of radius ``θ_r`` on a grid defined by ``[θ_x, θ_y]``. 
The center of the disk is at ``\\pmb{β} = (β_x, β_y)``. By default, the source profile is constant and every
pixel has a value of 1.0 and we can scale it using the amplitude ``A``. The corresponding formula is:
```math
S(θ_x, θ_y) =
\\begin{cases}
   A, & \\text{if } (θ_x - β_x)^2 + (θ_y - β_y)^2 ≤ θ_r^2 \\\\
   0, & \\text{otherwise}
\\end{cases}
```
"""
function disk(θx::Matrix{<:RV}, θy::Matrix{<:RV}, θr::RV, β::NTuple{2, RV}; A::RV=1.0)::Matrix{<:RV}
   # Initialize an empty source grid
   src = zero(θx)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - β[1]
         dy = θy[i, j] - β[2]
         src[i, j] = (dx^2 + dy^2 <= θr^2) ? A * 1.0 : 0.0
      end
   end
   return src
end


"""
    gaussian(θx::Matrix{<:RV}, θy::Matrix{<:RV}, σx::RV, σy::RV, β::NTuple{2, RV}; A::RV=1.0) --> Matrix{<:RV}

Creates a Gaussian source profile on a grid defined by ``[θ_x, θ_y]``. Standard deviations along 
``(x, y)`` axis are given by ``(σ_x, σ_y)``. The center of the Gaussian is at ``\\pmb{β} = (β_x, β_y)``. 
The overall normalization is determined by ``A``. The corresponding formula is:
```math
S(θ_x, θ_y) = \\frac{A}{2 π σ_x σ_y} \\exp\\left[-\\frac{1}{2} \\left(\\frac{(θ_x - β_x)^2}{σ_x^2} + \\frac{(θ_y - β_y)^2}{σ_y^2}\\right)\\right]
```
"""
function gaussian(θx::Matrix{<:RV}, θy::Matrix{<:RV}, σx::RV, σy::RV, β::NTuple{2, RV}; A::RV=1.0)::Matrix{<:RV}
   # Initialize an empty source grid
   src = zero(θx)

   # Normalization factor
   amplitude::Float64 = A / (2.0 * π * σx * σy)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - β[1]
         dy = θy[i, j] - β[2]
         src[i, j] = amplitude * exp(-0.5 * (dx^2 / σx^2 + dy^2 / σy^2))
      end
   end
   return src
end


@inline function b_n(n::RV)
   if 0.06 < n < 0.36
      return 0.01945 - n * (0.8902 - n * (10.95 - n * (19.67 - 13.47 * n)))
   else
      return 2.0*n - (1.0/3.0) + (4.0/405.0/n) + (46.0/25515/n^2) + (131.0/1148175.0/n^3) - (2194697.0/30690717750.0/n^4)
   end
end

"""
    sersic(θx::Matrix{<:RV}, θy::Matrix{<:RV}, n::RV, θe::RV, β::NTuple{2, RV}; A::RV=1.0) --> Matrix{<:RV}
Creates a Sersic source profile on a grid defined by ``[θ_x, θ_y]``. The Sersic index is given by ``n``
and the effective radius is given by ``θ_e``. The center of the Sersic profile is at ``\\pmb{β} = (β_x, β_y)``.
The overall normalization is determined by ``A``. The corresponding formula is:
```math
S(θ_x, θ_y) = \\frac{A \\,b_n^{2n}}{π θ_e^2 \\, Γ(2n+1)} 
\\exp\\left[-b_n \\left(\\frac{\\sqrt{(θ_x - β_x)^2 + (θ_y - β_y)^2}}{θ_e}\\right)^{1/n}\\right],
```
where,
```math
b_n = 
\\begin{cases}
0.01945 - 0.8902\\:n + 10.95\\:n^2 - 19.67\\:n^3 + 13.43\\:n^4, & 0.06 < n < 0.36 \\\\
2n - \\frac{1}{3} + \\frac{4}{405\\:n} + \\frac{46}{25515\\:n^2} + \\frac{131}{1148175\\:n^3} - \\frac{2194697}{30690717750\\:n^4}, & n > 0.36
\\end{cases}
```
"""
function sersic(θx::Matrix{<:RV}, θy::Matrix{<:RV}, n::RV, θe::RV, β::NTuple{2, RV}; A::RV=1.0)::Matrix{<:RV}
   # Initialize an empty source grid
   src = zero(θx)
   bn::Float64 = b_n(n)

   # Normalization factor
   amplitude::Float64 = A * bn^(2n) / (π * θe^2 * gamma(2.0 * n + 1.0))

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - β[1]
         dy = θy[i, j] - β[2]
         src[i, j] = amplitude * exp(-bn * (√(dx^2 + dy^2) / θe)^(1/n))
      end
   end
   return src
end

end