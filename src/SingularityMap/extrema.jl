function _extrema!(nex::Int64, rex::Matrix{Float64}, de1::Array{Float64, 3})
   # Local variables
   nx, ny = size(de1, 1), size(de1, 2)
   
   l   = 0
   k11 = 0
   k12 = 0
   k21 = 0
   k22 = 0

   xex1 = Matrix{Float64}(undef, Int(nx*ny/8), 2)
   yex1 = Matrix{Float64}(undef, Int(nx*ny/8), 2)
   xex2 = Matrix{Float64}(undef, Int(nx*ny/8), 2)
   yex2 = Matrix{Float64}(undef, Int(nx*ny/8), 2)

   is11 = zeros(Int64, Int(nx*ny/8))
   is12 = zeros(Int64, Int(nx*ny/8))
   is21 = zeros(Int64, Int(nx*ny/8))
   is22 = zeros(Int64, Int(nx*ny/8))


   ax1, ax2 = axes(ddf, 1)[4:end-3], axes(ddf, 2)[4:end-3]
   @inbounds for j in ax2
      @inbounds for i in ax1
         l  = 1
         j1 = i + 1
         if de1[j1, j, l] * de1[i, j, l] ≤ 0.0
            k11  = k11 + 1
            xex1[k11, 1] = float(i) - de1[i, j, l] / (de1[j1, j, l] - de1[i, j, l])
            yex1[k11, 1] = float(j)
         end
         
         j1 = j + 1
         if de1[i, j1, l] * de1[i, j, l] ≤ 0.0
            k12  = k12 + 1
            xex2[k12, 1] = float(i)
            yex2[k12, 1] = float(j) - de1[i, j, l] / (de1[i, j1, l] - de1[i, j, l])
         end
         
         l = 2
         j1 = i + 1
         if de1[j1, j, l] * de1[i, j, l] ≤ 0.0
            k21 = k21 + 1
            xex1[k21, 2] = float(i) - de1[i, j, l] / (de1[j1, j, l] - de1[i, j, l])
            yex1[k21, 2] = float(j)
         end
         
         j1 = j + 1
         if de1[i, j1, l] * de1[i, j, l] ≤ 0.0
            k22 = k22 + 1
            xex2[k22, 2] = float(i)
            yex2[k22, 2] = float(j) - de1[i, j, l] / (de1[i, j1, l] - de1[i, j, l])
         end
      end
   end

   nex = 0
   @inbounds for i in 1:k11
      @inbounds for j in 1:k22

      end
   end

   return nothing
end