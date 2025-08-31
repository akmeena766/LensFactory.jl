"""
    Sources
"""
module Sources


# LensFactory modules to use
using ..Constants

# Source profiles to export
export disk
export gaussian


"""
    disk(θx::Matrix{<:RV}, θy::Matrix{<:RV}, r::RV, β::NTuple{2, RV}; A::RV=1.0) --> Matrix{<:RV}

Creates a disk source profile of radius ``r`` on a grid defined by ``[θ_x, θ_y]``. 
The center of the disk is at ``\\pmb{β} = (β_x, β_y)``. By default, the source profile is constant and every
pixel has a value of 1.0 and we can scale it using the amplitude ``A``. The corresponding formula is:
```math
S(θ_x, θ_y) =
\\begin{cases}
   A, & \\text{if } (θ_x - β_x)^2 + (θ_y - β_y)^2 ≤ r^2 \\\\
   0, & \\text{otherwise}
\\end{cases}
```
"""
function disk(θx::Matrix{<:RV}, θy::Matrix{<:RV}, radius::RV, β::NTuple{2, RV}; A::RV=1.0)::Matrix{<:RV}
   # Initialize an empty source grid
   src = zero(θ_x)

   # Local variables for calculations
   dx::Float64 = 0.0
   dy::Float64 = 0.0

   ax1, ax2 = axes(θ_x, 1), axes(θ_x, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θ_x[i, j] - β[1]
         dy = θ_y[i, j] - β[2]
         src[i, j] = (dx^2 + dy^2 <= radius^2) ? A * 1.0 : 0.0
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

   # Local variables for calculations
   dx::Float64 = 0.0
   dy::Float64 = 0.0
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

end