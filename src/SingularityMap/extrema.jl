function _extrema!(nex::Int64, rex::Matrix{Float64}, de::Array{Float64, 3})
   # Local variables
   nx, ny = size(de, 2), size(de, 3)
   half_x = 0.5 * nx
   half_y = 0.5 * ny

   l::Int64   = 0
   k11::Int64 = 0
   k12::Int64 = 0
   k21::Int64 = 0
   k22::Int64 = 0

   is11 = zeros(Int64, round(Int64, nx*ny/8))
   is12 = zeros(Int64, round(Int64, nx*ny/8))
   is21 = zeros(Int64, round(Int64, nx*ny/8))
   is22 = zeros(Int64, round(Int64, nx*ny/8))

   xex1 = zeros(Float64, 2, round(Int64, nx*ny/8))
   yex1 = zeros(Float64, 2, round(Int64, nx*ny/8))
   xex2 = zeros(Float64, 2, round(Int64, nx*ny/8))
   yex2 = zeros(Float64, 2, round(Int64, nx*ny/8))

   ax1, ax2 = axes(de, 2)[4:end-3], axes(de, 3)[4:end-3]
   @inbounds for j in ax2
      @inbounds for i in ax1
         l  = 1
         i1 = i + 1
         if de[l, i1, j] * de[l, i, j] ≤ 0.0
            k11  = k11 + 1
            xex1[1, k11] = float(i) - de[l, i, j] / (de[l, i1, j] - de[l, i, j])
            yex1[1, k11] = float(j)
         end
         
         j1 = j + 1
         if de[l, i, j1] * de[l, i, j] ≤ 0.0
            k12  = k12 + 1
            xex2[1, k12] = float(i)
            yex2[1, k12] = float(j) - de[l, i, j] / (de[l, i, j1] - de[l, i, j])
         end
         
         l = 2
         i1 = i + 1
         if de[l, i1, j] * de[l, i, j] ≤ 0.0
            k21 = k21 + 1
            xex1[2, k21] = float(i) - de[l, i, j] / (de[l, i1, j] - de[l, i, j])
            yex1[2, k21] = float(j)
         end
         
         j1 = j + 1
         if de[l, i, j1] * de[l, i, j] ≤ 0.0
            k22 = k22 + 1
            xex2[2, k22] = float(i)
            yex2[2, k22] = float(j) - de[l, i, j] / (de[l, i, j1] - de[l, i, j])
         end
      end
   end

   @inbounds for i in 1:k11
      @inbounds for j in 1:k22
         if is22[j] == 0 && is11[i] == 0
            dx = abs(xex1[1, i] - xex2[2, j])
            if dx ≥ half_x
               dx = nx - dx
            end

            dy = abs(yex1[1, i] - yex2[2, j])
            if dy ≥ half_y
               dy = ny - dy
            end

            if dx < 1.0 && dy < 1.0
               nex = nex + 1
               rex[1, nex] = xex1[1, i]
               rex[2, nex] = yex1[1, i]
               is22[j] = 1
               is11[i] = 1
            end
         end
      end
   end

   @inbounds for i in 1:k12
      @inbounds for j in 1:k21
         if is21[j] == 0 && is12[i] == 0
            dx = abs(xex2[1, i] - xex1[2, j])
            if dx ≥ half_x
               dx = nx - dx
            end

            dy = abs(yex2[1, i] - yex1[2, j])
            if dy ≥ half_y
               dy = ny - dy
            end

            if dx < 1.0 && dy < 1.0
               nex = nex + 1
               rex[1, nex] = xex2[1, i]
               rex[2, nex] = yex2[1, i]
               is21[j] = 1
               is12[i] = 1
            end
         end
      end
   end
   return nex
end


function _maxima!(nmax::Int64, rmax::T, e1::T, nex::Int64, rex::T) where T <: Matrix{Float64}
   nx, ny = size(e1)
   half_x = 0.5 * nx
   half_y = 0.5 * ny

   ntemp = 0
   rtemp = zeros(Float64, 2, round(Int64, nx*ny/8))
   @inbounds for ii in 1:nex
      i0 = round(Int64, rex[1, ii])
      j0 = round(Int64, rex[2, ii])
      for k1 in i0-1:i0+1
         for k2 in j0-1:j0+1
            i = k1
            j = k2
            l = 0
            for i1 in -2:2
               for i2 in -2:2
                  if i1 ≠ 0 || i2 ≠ 0
                     j1 = i1 + i
                     j2 = i2 + j
                     if e1[i, j] > e1[j1, j2]
                        l = l + 1
                     end
                  end
               end
            end
            if l == 24
               ntemp = ntemp + 1
               rtemp[1, ntemp] = float(i)
               rtemp[2, ntemp] = float(j)
            end
         end
      end
   end

   i = 1
   nmax = 1
   for l in 1:2
      rmax[l, nmax] = rtemp[l, i]
   end
   for i in 2:ntemp
      ii = 0
      for j in 1:i-1
         d = 0.0
         for l in 1:2
            dx = abs(rtemp[l, i] - rtemp[l, j])
            if dx ≥ half_x
               dx = nx - dx
            end
            d = d + dx^2
         end
         if d < 2.0
            ii = 1
         end
      end
      if ii == 0
         nmax = nmax + 1
         for l in 1:2
            rmax[l, nmax] = rtemp[l, i]
         end
      end
   end
   return nmax
end


function _minima!(nmin::Int64, rmin::T, e1::T, nex::Int64, rex::T) where T <: Matrix{Float64}
   nx, ny = size(e1)
   half_x = 0.5 * nx
   half_y = 0.5 * ny

   ntemp = 0
   rtemp = zeros(Float64, 2, round(Int64, nx*ny/8))
   @inbounds for ii in 1:nex
      i0 = round(Int64, rex[1, ii])
      j0 = round(Int64, rex[2, ii])
      for k1 in i0-1:i0+1
         for k2 in j0-1:j0+1
            i = k1
            j = k2
            l = 0
            for i1 in -2:2
               for i2 in -2:2
                  if i1 ≠ 0 || i2 ≠ 0
                     j1 = i1 + i
                     j2 = i2 + j
                     if e1[i, j] > e1[j1, j2]
                        l = l + 1
                     end
                  end
               end
            end
            if l == 24
               ntemp = ntemp + 1
               rtemp[1, ntemp] = float(i)
               rtemp[2, ntemp] = float(j)
            end
         end
      end
   end

   i = 1
   nmin = 1
   for l in 1:2
      rmin[l, nmin] = rtemp[l, i]
   end
   for i in 2:ntemp
      ii = 0
      for j in 1:i-1
         d = 0.0
         for l in 1:2
            dx = abs(rtemp[l, i] - rtemp[l, j])
            if dx ≥ half_x
               dx = nx - dx
            end
            d = d + dx^2
         end
         if d < 2.0
            ii = 1
         end
      end
      if ii == 0
         nmin = nmin + 1
         for l in 1:2
            rmin[l, nmin] = rtemp[l, i]
         end
      end
   end
   return nmin
end