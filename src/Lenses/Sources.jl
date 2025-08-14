"""
    Sources
Testing :-)
"""
module Sources


# LensFactory modules to use
using ..Constants

# Source profiles to export
export disk
export gaussian


"""
    disk(θx::ROA, θy::ROA, r::Float64, β::NTuple{2, RV}; A::RV=1.0) --> Matrix{<:RV}

Creates a disk source profile of radius ``r`` on a grid defined by ``[θ_x, θ_y]``. 
``β`` is the center of the disk. By default, the source profile is constant and every
pixel has a value of 1.0 and we can scale it using the amplitude ``A``.
"""
function disk(θx::ROA, θy::ROA, radius::Float64, β::NTuple{2, RV}; A::RV=1.0)::Matrix{<:RV}
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
    gaussian(θx::ROA, θy::ROA, σx::RV, σy::RV, β::NTuple{2, RV}; A::RV=1.0) --> Matrix{<:RV}

Creates a Gaussian source profile on a grid defined by ``[θ_x, θ_y]``. Standard deviations along 
``(x, y)`` axis are given by ``(σ_x, σ_y)``. The center of the Gaussian is at ``β``. The overall 
normalization is determined by ``A``.
"""
function gaussian(θ_x::ROA, θ_y::ROA, σ_x::RV, σ_y::RV, β::NTuple{2, RV}; A::RV=1.0)::Matrix{<:RV}
   # Initialize an empty source grid
   src = zero(θ_x)

   # Local variables for calculations
   dx::Float64 = 0.0
   dy::Float64 = 0.0
   amplitude::Float64 = A / (2.0 * π * σ_x * σ_y)

   ax1, ax2 = axes(θ_x, 1), axes(θ_x, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θ_x[i, j] - β[1]
         dy = θ_y[i, j] - β[2]
         src[i, j] = amplitude * exp(-0.5 * (dx^2 / σ_x^2 + dy^2 / σ_y^2))
      end
   end
   return src
end

end