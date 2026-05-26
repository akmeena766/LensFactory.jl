function _a4points!(na4::S, ra4::T, e::T, na3::S, ra3::T, nex::S, rex::T, nmax::S, rmax::T, nmin::S, rmin::T) where {T <: Matrix{Float64}, S <: Int64}
   # Local variables
   nx, ny = size(e, 1), size(e, 2)
   half_x = 0.5 * nx
   half_y = 0.5 * ny

   nd4::Int64 = 0
   nt::Int64  = 0
   l1::Int64  = 0
   l2::Int64  = 0
   ie1::Int64 = 0
   ie2::Int64 = 0
   isv  = zeros(Int64,   round(Int64, nx*ny/8))
   ev   = zeros(Float64, round(Int64, nx*ny/8))
   rd4  = zeros(Float64, 2, round(Int64, nx*ny/8))
   link = zeros(Int64,   2, round(Int64, nx*ny/8))
   rt   = zeros(Float64, 2, round(Int64, nx*ny/8))

   for i in 1:na3
      j1 = floor(Int64, ra3[1, i])
      j2 = floor(Int64, ra3[2, i])
      dx = ra3[1, i] - float(j1)
      dy = ra3[2, i] - float(j2)
      
      ev[i] = (1.0 - dx - dy) * e[j1, j2] + dx * e[j1+1, j2] + dy * e[j1, j2+1]

      l1 = 0
      l2 = 0
      if dx < 1.0E-5
         for j in 1:na3
            if j ≠ i
               if ra3[2, j] ≥ float(j2) && ra3[2, j] ≤ float(j2+1)
                  if ra3[1, j] ≥ float(j1-1) && ra3[1, j] ≤ float(j1)
                     l1  = l1 +1
                     ie1 = j
                  end
               end

               if ra3[2, j] ≥ float(j2) && ra3[2, j] ≤ float(j2+1)
                  if ra3[1, j] ≥ float(j1) && ra3[1, j] ≤ float(j1+1)
                     l2  = l2 +1
                     ie2 = j
                  end
               end
            end
         end
      elseif dy < 1.0E-5
         for j in 1:na3
            if j ≠ i
               if ra3[1, j] ≥ float(j1) && ra3[1, j] ≤ float(j1+1)
                  if ra3[2, j] ≥ float(j2-1) && ra3[2, j] ≤ float(j2)
                     l1  = l1 + 1
                     ie1 = j
                  end
               end

               if ra3[1, j] ≥ float(j1) && ra3[1, j] ≤ float(j1+1)
                  if ra3[2, j] ≥ float(j2) && ra3[2, j] ≤ float(j2+1)
                     l2  = l2 + 1
                     ie2 = j
                  end
               end
            end
         end
      end

      k = l1 + l2
      isv[i] = k
      if l1 == 1
         link[1, i] = ie1
      end
      if l2 == 1
         link[2, i] = ie2
      end
      if isv[i] > 2
         nd4 = nd4 + 1
      end
   end

   nd4 = 0
   for i in 1:na3
      if isv[i] == 2
         j1 = link[1, i]
         j2 = link[2, i]
         if j1 ≠ 0 && j2 ≠ 0
            k1 = link[1, j1]
            if i == k1
               k1 = link[2, j1]
            end
            k2 = link[1, j2]
            if i == k2
               k2 = link[2, j2]
            end
            if k1 ≠ 0 && k2 ≠ 0
               l1 = link[1, k1]
               if l1 == j1
                  l1 = link[2, k1]
               end
               l2 = link[1, k2]
               if l2 == j2
                  l2 = link[2, k2]
               end
               if l1 ≠ 0 && l2 ≠ 0
                  if ev[i] ≥ ev[j1] && ev[i] ≥ ev[j2]
                     if ev[i] ≥ ev[k1] && ev[i] ≥ ev[k2]
                        if ev[i] ≥ ev[l1] && ev[i] ≥ ev[l2]
                           ia4 = 0
                           for k in 1:nmax
                              d = 0.0
                              
                              dx = abs(rmax[1, k] - ra3[1, i])
                              if dx > half_x
                                 dx = nx - dx
                              end
                              d = d + dx * dx

                              dy = abs(rmax[2, k] - ra3[2, i])
                              if dy > half_y
                                 dy = ny - dy
                              end
                              d = d + dy * dy

                              if d ≤ 2.0
                                 ia4 = 1
                              end
                           end

                           if ia4 == 0
                              for k in 1:nex
                                 d = 0.0
                                 
                                 dx = abs(rex[1, k] - ra3[1, i])
                                 if dx > half_x
                                    dx = nx - dx
                                 end
                                 d = d + dx * dx

                                 dy = abs(rex[2, k] - ra3[2, i])
                                 if dy > half_y
                                    dy = ny - dy
                                 end
                                 d = d + dy * dy

                                 if d ≤ 2.0
                                    nd4 = nd4 + 1
                                    rd4[1, nd4] = ra3[1, i]
                                    rd4[2, nd4] = ra3[2, i]
                                    ia4 = 1
                                 end
                              end
                           end

                           if ia4 == 0
                              na4 = na4 + 1
                              ra4[1, na4] = ra3[1, i]
                              ra4[2, na4] = ra3[2, i]
                           end
                        end
                     end
                  end

                  if ev[i] ≤ ev[j1] && ev[i] ≤ ev[j2]
                     if ev[i] ≤ ev[k1] && ev[i] ≤ ev[k2]
                        if ev[i] ≤ ev[l1] && ev[i] ≤ ev[l2]
                           ia4 = 0
                           for k in 1:nmin
                              d = 0.0

                              dx = abs(rmin[1, k] - ra3[1, i])
                              if dx > half_x
                                 dx = nx - dx
                              end
                              d = d + dx * dx

                              dy = abs(rmin[2, k] - ra3[2, i])
                              if dy > half_y
                                 dy = ny - dy
                              end
                              d = d + dy * dy
                              
                              if d ≤ 2.0
                                 ia4 = 1
                              end
                           end

                           if ia4 == 0
                              for k in 1:nex
                                 d = 0.0
                                 
                                 dx = abs(rex[1, k] - ra3[1, i])
                                 if dx > half_x
                                    dx = nx - dx
                                 end
                                 d = d + dx * dx

                                 dy = abs(rex[2, k] - ra3[2, i])
                                 if dy > half_y
                                    dy = ny - dy
                                 end
                                 d = d + dy * dy

                                 if d ≤ 2.0
                                    nt = nt + 1
                                    rt[1, nt] = ra3[1, i]
                                    rt[2, nt] = ra3[2, i]
                                 end
                              end
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end
   return na4
end