function _a3lines!(na3::Int64, ra3::T, e::T, de::S, q::S) where {T <: Matrix{Float64}, S <: Array{Float64,3}}
   # Local variables
   nx, ny = size(e)

   dot1::Float64 = 0.0
   dot2::Float64 = 0.0
   sign::Float64 = 0.0

   ax1, ax2 = axes(e, 1)[10:end-9], axes(e, 2)[10:end-9]
   @inbounds for j in ax2
      @inbounds for i in ax1
         k = i + 1
         dot1 = 0.0
         for l in 1:2
            dot1 = dot1 + q[l, i, j] * q[l, k, j]
         end
         if dot1 < 0.0
            sign = -1.0
         else
            sign = +1.0
         end
         dot1 = 0.0
         dot2 = 0.0
         for l in 1:2
            dot1 = dot1 +        q[l, i, j] * de[l, i, j]
            dot2 = dot2 + sign * q[l, k, j] * de[l, k, j]
         end
         if dot1 * dot2 ≤ 0.0
            na3 = na3 + 1
            ra3[1, na3] = float(i) - dot1 / (dot2 - dot1)
            ra3[2, na3] = float(j)
         end

         k = j + 1
         dot1 = 0.0
         for l in 1:2
            dot1 = dot1 + q[l, i, j] * q[l, i, k]
         end
         if dot1 < 0.0
            sign = -1.0
         else
            sign = +1.0
         end
         dot1 = 0.0
         dot2 = 0.0
         for l in 1:2
            dot1 = dot1 +        q[l, i, j] * de[l, i, j]
            dot2 = dot2 + sign * q[l, i, k] * de[l, i, k]
         end
         if dot1 * dot2 ≤ 0.0
            na3 = na3 + 1
            ra3[1, na3] = float(i)
            ra3[2, na3] = float(j) - dot1 / (dot2 - dot1)
         end
      end
   end
   return na3
end