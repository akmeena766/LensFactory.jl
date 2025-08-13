module Sources


# LensFactory modules to use
using ..Constants

# Source profiles to export
export gaussian

function gaussian(θ_x::ROA, θ_y::ROA, σ_x::RV, σ_y::RV, β::NTuple{2, RV}; A::RV=1.0)
    # Initialize an empty source grid
    src = zero(θ_x)
    dx::Float64 = 0.0
    dy::Float64 = 0.0
    amplitude::Float64 = A / (2π * σ_x * σ_y)

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