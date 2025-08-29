module PolygonOps

export shoelace
export hao_sun

function shoelace(curve::Vector{Vector{Float64}})::Float64
   # Close the polygon if needed
   if curve[1] != curve[end]
      push!(curve, curve[1])
   end

   # Calculate the area
   area = 0.0
   for i in 1:(length(curve)-1)
      x1, y1 = curve[i]
      x2, y2 = curve[i+1]
      area += x1*y2 - x2*y1
   end
   return 0.5 * abs(area)
end


function hao_sun(point, polygon)
   """
   Algorithm to find determine if a point is in the given polygon.
   Taken from Hao and Sun (2018): https://doi.org/10.3390/sym10100477
   outputs:
      -1: on the edge
       0: outside
       1: inside
   """
   
   k = 0

   xp = point[1]
   yp = point[2]

   @inbounds for i in 1:length(polygon)-1
      v1 = polygon[i][2] - yp
      v2 = polygon[i+1][2] - yp

      if v1 < 0 && v2 < 0 || v1 > 0 && v2 > 0
         continue
      end

      u1 = polygon[i][1] - xp
      u2 = polygon[i+1][1] - xp

      f = (u1 * v2) - (u2 * v1)

      if v2 > 0 && v1 <= 0
         if f > 0
            k += 1
         elseif f == 0
            return -1
         end
      elseif v1 > 0 && v2 <= 0
         if f < 0
            k += 1
         elseif f == 0
            return -1
         end
      elseif v2 == 0 && v1 < 0
         if f == 0
            return -1
         end   
      elseif v1 == 0 && v2 < 0
         if f == 0
            return -1
         end 
      elseif v1 == 0 && v2 == 0
         if u2 <= 0 && u1 >= 0
            return on
         elseif u1 <= 0 && u2 >= 0
            return on
         end
      end
   end
   
   if k % 2 == 0
      return 0
   end

   return 1
end

end