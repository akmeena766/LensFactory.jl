function _extrema!(nex::Int64, rex::Matrix{Float64}, de1::Array{Float64, 3})
   # Local variables
   nx, ny = size(de1, 1), size(de1, 2)
   half_x = 0.5 * nx
   half_y = 0.5 * ny
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

   nex = 0
   rex = Matrix{Float64}(undef, Int(nx*ny/8), 2)

   ax1, ax2 = axes(de1, 1)[4:end-3], axes(de1, 2)[4:end-3]
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

   @inbounds for i in 1:k11
      @inbounds for j in 1:k22
         if is22[j] == 0 && is11[i] == 0
            dx = abs(xex1[i, 1] - xex2[j, 2])
            if dx ≥ half_x
               dx = nx - dx
            end

            dy = abs(yex1[i, 1] - yex2[j, 2])
            if dy ≥ half_y
               dy = ny - dy
            end

            if dx < 1.0 && dy < 1.0
               nex = nex + 1
               rex[nex, 1] = xex1[i, 1]
               rex[nex, 2] = yex1[i, 1]
               is22[j] = 1
               is11[i] = 1
            end
         end
      end
   end

   @inbounds for i in 1:k12
      @inbounds for j in 1:k21
         if is21[j] == 0 && is12[i] == 0
            dx = abs(xex2[i, 1] - xex1[j, 2])
            if dx ≥ half_x
               dx = nx - dx
            end

            dy = abs(yex2[i, 1] - yex1[j, 2])
            if dy ≥ half_y
               dy = ny - dy
            end

            if dx < 1.0 && dy < 1.0
               nex = nex + 1
               rex[nex, 1] = xex2[i, 1]
               rex[nex, 2] = yex2[i, 1]
               is21[j] = 1
               is12[i] = 1
            end
         end
      end
   end
   return nothing
end


function _minima!()
   
end


function _maxima!()
   
end