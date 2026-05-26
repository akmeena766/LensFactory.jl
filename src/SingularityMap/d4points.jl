function _d4points(nd4::Int64, rd4::Matrix{Float64}, ddf::Array{Float64,3})
   # Local variables
   nx, ny = size(ddf, 2), size(ddf, 3)

   ntemp::Int64 = 0
   n1122::Int64 = 0
   n1212::Int64 = 0

   x1 = zeros(Float64, 4)
   x2 = zeros(Float64, 4)
   y1 = zeros(Float64, 4)
   y2 = zeros(Float64, 4)
   rtemp = zeros(Float64, 2, round(Int64, nx*ny/8))

   ax1, ax2 = axes(ddf, 2)[10:end-9], axes(ddf, 3)[10:end-9]
   @inbounds for j in ax2
      @inbounds for i in ax1
         n1122 = 0
         n1212 = 0

         k1 = i + 1
         k2 = j + 1

         diff1 = ddf[1, i,   j] - ddf[2, i,   j]
         diff2 = ddf[1, k1,  j] - ddf[2, k1,  j]
         diff3 = ddf[1, i,  k2] - ddf[2, i,  k2]
         diff4 = ddf[1, k1, k2] - ddf[2, k1, k2]
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

         if ddf[3, i, j] * ddf[3, k1, j] ≤ 0.0
            if abs(ddf[3, i, j] - ddf[3, k1, j]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i) + ddf[3, i, j] / (ddf[3, i, j] - ddf[3, k1, j])
               y2[n1212] = float(j)
            end
         end
         if ddf[3, i, k2] * ddf[3, k1, k2] ≤ 0.0
            if abs(ddf[3, i, k2] - ddf[3, k1, k2]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i) + ddf[3, i, k2] / (ddf[3, i, k2] - ddf[3, k1, k2])
               y2[n1212] = float(j + 1)
            end
         end
         if ddf[3, i, j] * ddf[3, i, k2] ≤ 0.0
            if abs(ddf[3, i, j] - ddf[3, i, k2]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i)
               y2[n1212] = float(j) + ddf[3, i, j] / (ddf[3, i, j] - ddf[3, i, k2])
            end
         end
         if ddf[3, k1, j] * ddf[3, k1, k2] ≤ 0.0
            if abs(ddf[3, k1, j] - ddf[3, k1, k2]) > 1.0E-15
               n1212 = n1212 + 1
               y1[n1212] = float(i + 1)
               y2[n1212] = float(j) + ddf[3, k1, j] / (ddf[3, k1, j] - ddf[3, k1, k2])
            end
         end

         if n1212 == 2 && n1122 == 2
            ntemp = ntemp + 1

            # Direction vectors for n1122 contour
            dx1 = x1[2] - x1[1]
            dx2 = x2[2] - x2[1]

            # Direction vectors for n1212 contour
            dy1 = y1[2] - y1[1]
            dy2 = y2[2] - y2[1]

            # Solve: (x1[1], x2[1]) + t*(dx1, dx2) = (y1[1], y2[1]) + s*(dy1, dy2)
            denom = dx1 * dy2 - dx2 * dy1

            if abs(denom) > 1.0E-15
               # contours not parallel
               t = ((y1[1] - x1[1]) * dy2 - (y2[1] - x2[1]) * dy1) / denom
               rtemp[1, ntemp] = x1[1] + t * dx1
               rtemp[2, ntemp] = x2[1] + t * dx2
            else
               # contours parallel (fallback to centroid)
               rtemp[1, ntemp] = (x1[1] + x1[2] + y1[1] + y1[2]) / 4.0
               rtemp[2, ntemp] = (x2[1] + x2[2] + y2[1] + y2[2]) / 4.0
            end
         end
      end
   end

   # Remove duplicates
   @inbounds for i in 1:ntemp
      k1 = 0
      @inbounds for j in i+1:ntemp
         dr = 0.0
         for l in 1:2
            dr = dr + (rtemp[l, j] - rtemp[l, i])^2
         end
         if dr < 1.0
            k1 = 1
         end
      end
      if k1 == 0
         nd4 = nd4 + 1
         for l in 1:2
            rd4[l, nd4] = rtemp[l, i]
         end
      end
   end
   return nd4
end