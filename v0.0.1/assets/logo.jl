using Luxor
using Colors

function micro_mu()
	# Microlensing magnification function
	A(u) = (u^2 + 2) / (u * sqrt(u^2 + 4))

	# Parameters
	u0 = 0.1        # impact parameter
	tE = 1.0        # Einstein timescale
	t0 = 0.0        # time of peak

	# Time range
	ts = range(-2, 2, length=100)   # -2 tE to +2 tE
	us = sqrt.(u0^2 .+ ((ts .- t0) ./ tE).^2)
	amps = A.(us)
	return ts, amps
end


function strong_image()
	# Define grid
   πxel_size = 0.002
	s = collect(-5:πxel_size:5)
	xg = zeros(length(s), length(s))
	yg = zeros(length(s), length(s))

	src = zeros(length(s), length(s))
	img = zeros(length(s), length(s))

	for j=1:size(s,1)
    	for i=1:size(s,1)
        	xg[i, j] = s[i]
        	yg[i, j] = s[j]
         if sqrt((xg[i, j]-0.1)^2 + (yg[i, j]-0.1)^2) <= 0.1
			   src[i, j] = 1.0
         end
    	end
	end

   for j=1:length(s)
    	for i=1:length(s)
         xg[i, j] = xg[i, j] + 1E-20
         yg[i, j] = yg[i, j] + 1E-20

         y1 = xg[i, j] - xg[i, j] / √( xg[i, j]^2 + yg[i, j]^2 )
         y2 = yg[i, j] - yg[i, j] / √( xg[i, j]^2 + yg[i, j]^2 )

         πxel_x = round(Int64, y1/πxel_size + 0.5 * length(s) + 1.0)
         πxel_y = round(Int64, y2/πxel_size + 0.5 * length(s) + 1.0)

         if (1 <= πxel_x <= length(s)) && (1 <= πxel_y <= length(s))
            img[i, j] = src[πxel_x, πxel_y]
         end
      end
   end

   return xg, yg, img
end


function create_lensfactory_logo()
	R = 100
	Drawing(7R/2, 7R/2, "logo.png")
	origin()
	
	translate(0, R/4)
	juliacircles(R)
    
   # Center of the drawing and circle parameters
   center_point = Point(7R/4, 7R/4)
   
   # Calculate positions for three circles in Julia's triangular arrangement
   position1 = Point(0, -4R/4)
	
	xv = position1[1] * cos(120*π/180) - position1[2] * sin(120*π/180)
	yv = position1[1] * sin(120*π/180) + position1[2] * cos(120*π/180)
	position2 = Point(xv, yv)
	
	xv = position1[1] * cos(120*π/180) + position1[2] * sin(120*π/180)
	yv = position1[1] * sin(120*π/180) + position1[2] * cos(120*π/180)
	position3 = Point(xv, yv)


   # setcolor(0.22, 0.596, 0.149, 1.0)
	# circle(position1, 3R/4, :stroke)

   # sethue(RGB(0.584, 0.345, 0.698))
	# circle(position2, 3R/4, :stroke)
	
   # sethue(RGB(0.796, 0.235, 0.2))
   # circle(position3, 3R/4, :stroke)

	# Strong lensing -- upper circle
	xg, yg, img = strong_image()

   function scale_point1(x, y)
      x = +x*54 + position1[1]
	   y = -y*54 + position1[2] #+ 50
	   return Point(x, y)
   end

   sethue(RGB(1.0, 1.0, 1.0))#sethue(RGB(0.22, 0.596, 0.149))
   for j in axes(img,2)
      for i in axes(img,1)
         if img[i, j] == 1 
            circle(scale_point1(xg[i, j], yg[i, j]), 1, :fill)
         end
      end
   end

   # sp1 = Point(position1[1]-3R/4, position1[2])
   # ep1 = Point(position1[1]+3R/4, position1[2])
   
   # setdash("dashed")
   # line(sp1, ep1, :stroke)

   # angle = -90
   # xv = position1[1] - 3R/4 * cos(angle*π/180)
	# yv = position1[2] - 3R/4 * sin(angle*π/180)
   # sp2 = Point(xv, yv)

   # xv = position1[1] + 3R/4 * cos(angle*π/180)
	# yv = position1[2] + 3R/4 * sin(angle*π/180)
   # ep2 = Point(xv, yv)

   # line(sp2, ep2, :stroke)

	# Microlensing -- lower right circle
	# Draw microlensing light curve
	ts, amps = micro_mu()
	
	function scale_point2(x, y)
	   x = +x*28 + position2[1]
	   y = -y*11 + position2[2] + 50
	   return Point(x, y)
	end
	
   setdash("solid")
	sethue(RGB(1.0, 1.0, 1.0))#sethue(RGB(0.584, 0.345, 0.698))
	setline(5)
	
	points = [scale_point2(t, y) for (t, y) in zip(ts, amps)]
	poly(points, :stroke)

   # Weak lensing
   sethue(RGB(1.0, 1.0, 1.0))#sethue(RGB(0.796, 0.235, 0.2))
   setline(2)
   # Draw ellipses at regular intervals around the circle
    num_ellipses = 12
    angles = range(0, 360, length=num_ellipses+1)[1:end-1]  # 0 to 330 degrees
    
    for (i, angle) in enumerate(angles)
      # Calculate position on circle
      pos = position3 + polar(R/2, deg2rad(angle))
      
      @layer begin
         # Move to ellipse position and rotate to face outward
         translate(pos)
         rotate(deg2rad(angle/2))
         
         # Draw ellipse with alternating colors
         ellipse(Point(0, 0), 20, 10, :stroke)
      end
    end

   finish()
   preview()
end

# Create the logo
create_lensfactory_logo()