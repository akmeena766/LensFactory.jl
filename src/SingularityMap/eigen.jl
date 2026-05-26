function _eigen!(e1::T, e2::T, q1::S, q2::S, ddf::S) where {T <: Matrix{Float64}, S <: Array{Float64,3}}
   # Initilize local variables
   tidal = zeros(Float64, 2, 2)
   
   # Calculate eigenvalues and eigenvectors
   ax1, ax2 = axes(ddf, 2)[4:end-3], axes(ddf, 3)[4:end-3]
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Get the jacobian matrix components
         tidal[1, 1] = ddf[1, i, j]
         tidal[2, 2] = ddf[2, i, j]
         tidal[1, 2] = ddf[3, i, j]
         tidal[2, 1] = ddf[3, i, j]

         term1 = 0.5 * (tidal[1, 1] + tidal[2, 2])
         term2 = 0.5 * sqrt( (tidal[1, 1] + tidal[2, 2])^2 -  4.0*(tidal[1, 1]*tidal[2, 2] - tidal[1, 2]*tidal[2, 1]) )

         a1 = term1 + term2
         a2 = term1 - term2

         # Update eigenvalues
         e1[i, j] = max(a1, a2)
         e2[i, j] = min(a1, a2)

         # Eigenvector for e1
         term1 = tidal[1, 1] - e1[i, j]
         term2 = tidal[2, 2] - e1[i, j]
         if term1 != 0.0
            e1x = -tidal[1, 2] / term1
            q1[1, i, j] = +e1x / sqrt(1.0 + e1x^2)
            q1[2, i, j] = +1.0 / sqrt(1.0 + e1x^2)
         else
            e1x = -tidal[2, 1] / term2
            q1[1, i, j] = +1.0 / sqrt(1.0 + e1x^2)
            q1[2, i, j] = +e1x / sqrt(1.0 + e1x^2)
         end

         # Eigenvector for e2
         term1 = tidal[1, 1] - e2[i, j]
         term2 = tidal[2, 2] - e2[i, j]
         if term1 != 0.0
            e2x = -tidal[1, 2] / term1
            q2[1, i, j] = +e2x / sqrt(1.0 + e2x^2)
            q2[2, i, j] = +1.0 / sqrt(1.0 + e2x^2)
         else
            e2x = -tidal[2, 1] / term2
            q2[1, i, j] = +1.0 / sqrt(1.0 + e2x^2)
            q2[2, i, j] = +e2x / sqrt(1.0 + e2x^2)
         end
      end
   end
   return nothing
end

function _gradeigen!(de1::S, de2::S, e1::T, e2::T) where {T <: Matrix{Float64}, S <: Array{Float64, 3}}
   ax1, ax2 = axes(e1, 1)[4:end-3], axes(e1, 2)[4:end-3]
   @inbounds for j in ax2
      @inbounds for i in ax1
         de1[1, i, j] = 0.5 * (e1[i+1, j] - e1[i-1, j])
         de1[2, i, j] = 0.5 * (e1[i, j+1] - e1[i, j-1])
         de2[1, i, j] = 0.5 * (e2[i+1, j] - e2[i-1, j])
         de2[2, i, j] = 0.5 * (e2[i, j+1] - e2[i, j-1])
      end
   end

   # Normalize the eigen values
   de1mag::Float64 = 0.0
   de2mag::Float64 = 0.0
   @inbounds for j in ax2
      @inbounds for i in ax1
         de1mag = 0.0
         de2mag = 0.0
         for l in 1:2
            de1mag = de1mag + de1[l, i, j] * de1[l, i, j]
            de2mag = de2mag + de2[l, i, j] * de2[l, i, j]
         end   
         de1mag=sqrt(de1mag)
         de2mag=sqrt(de2mag)
         if de1mag ≠ 0.0 && de2mag ≠ 0.0
            for l in 1:2
               de1[l, i, j] = de1[l, i, j] / de1mag
               de2[l, i, j] = de2[l, i, j] / de2mag
            end
         end
      end
   end
   return nothing
end