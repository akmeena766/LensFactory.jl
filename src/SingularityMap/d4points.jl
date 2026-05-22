function d4points(nd4::Int64, rd4::Matrix{Float64}, ddf::Array{Float64, 3})
   # Local variables
   nx, ny = size(e, 1), size(e, 2)
   
   ntemp::Int64 = 0
   n1122::Int64 = 0
   n1212::Int64 = 0
   
   x1 = zeros(Float64, 4)
   x2 = zeros(Float64, 4)
   y1 = zeros(Float64, 4)
   y2 = zeros(Float64, 4)
   rtemp = zeros(Float64, round(Int64, nx*ny/8), 2)

   ax1, ax2 = axes(e, 1)[10:end-9], axes(e, 2)[10:end-9]
   @inbounds for j in ax2
      @inbounds for i in ax1
         n1122 = 0
         n1212 = 0

         k1 = i + 1
         k2 = j + 1

         diff1 = ddf[i,   j, 1] - ddf[i,   j, 2]
         diff2 = ddf[k1,  j, 1] - ddf[k1,  j, 2]
         diff3 = ddf[i,  k2, 1] - ddf[i,  k2, 2]
         diff4 = ddf[k1, k2, 1] - ddf[k1, k2, 2]
         if diff1 * diff2 ≤ 0.0
            if abs(diff1 - diff2) > 1.0E-15
               n1122 = n1122 + 1
               x1[n1122] = float(i) + diff1 / (diff1 - diff2)
               x2[n1122] = float(j)
            end
         end
         if diff1 * diff3 ≤ 0.0
            if abs(diff1 - diff3) > 1.0E-15
               n1122 = n1122 + 1
               x1[n1122] = float(i)
               x2[n1122] = float(j) + diff1 / (diff1 - diff3)
            end
         end
         if diff2 * diff4 ≤ 0.0
            if abs(diff2 - diff4) > 1.0E-15
               n1122 = n1122 + 1
               x1[n1122] = float(i + 1)
               x2[n1122] = float(j) + diff2 / (diff2 - diff4)
            end
         end
         if diff3 * diff4 ≤ 0.0
            if abs(diff3 - diff4) > 1.0E-15
               n1122 = n1122 + 1
               x1[n1122] = float(i) + diff3 / (diff3 - diff4)
               x2[n1122] = float(j + 1)
            end
         end

         if ddf[i, j, 3] * ddf[k1, j, 3] ≤ 0.0
            if abs(ddf[i, j, 3] - ddf[k1, j, 3]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i) + ddf[i, j, 3] / (ddf[i, j, 3] - ddf[k1, j, 3])
               y2[n1212] = float(j)
            end
         end
         if ddf[i, k2, 3] * ddf[k1, k2, 3] ≤ 0.0
            if abs(ddf[i, k2, 3] - ddf[k1, k2, 3]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i) + ddf[i, k2, 3] / (ddf[i, k2, 3] - ddf[k1, k2, 3])
               y2[n1212] = float(j + 1)
            end
         end
         if ddf[i, j, 3] * ddf[i, k2, 3] ≤ 0.0
            if abs(ddf[i, j, 3] - ddf[i, k2, 3]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i)
               y2[n1212] = float(j) + ddf[i, j, 3] / (ddf[i, j, 3] - ddf[i, k2, 3])
            end
         end
         if ddf[k1, j, 3] * ddf[k1, k2, 3] ≤ 0.0
            if abs(ddf[k1, j, 3] - ddf[k1, k2, 3]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i + 1)
               y2[n1212] = float(j) + ddf[k1, j, 3] / (ddf[k1, j, 3] - ddf[k1, k2, 3])
            end
         end

         if n1212 == 2 && n1122 == 2
            ntemp = ntemp + 1
            rtemp(ntemp, 1) = float(i)
            rtemp(ntemp, 2) = float(j)
         end
      end
   end

   @inbounds for i in 1:ntemp
      k1 = 0
      @inbounds for j in i+1:ntemp
         dr = 0.0
         for l in 1:2
            dr = dr + (rtemp[j, l] - rtemp[i, l])^2
         end
         if dr < 4.0
            k1 = 1
         end
      end
      if k1 == 0
         nd4 = nd4 + 1
         for l in 1:2
            rd4[nd4, l] = rtemp[i, l]
         end
      end
   end
   return nd4
end