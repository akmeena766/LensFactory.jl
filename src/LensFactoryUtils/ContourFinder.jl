module ContourFinder


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# Various function to export
# --------------------------------------------------------------------------------------------------
export get_contour


function get_contour(x::AbstractMatrix{<:Real}, y::AbstractMatrix{<:Real}, z::AbstractMatrix{<:Real}, level::Real)
   # Checking the type of input data
   if !(axes(x) == axes(y) == axes(z))
      throw(ArgumentError("Incompatible input axes for x, y, and z."))
   end

   # Promoting all input data to same type
   ET = promote_type(map(eltype, (x, y, z))...)
   ET = ET <: Integer ? Float64 : ET

   # Get the cells with given level
   level_cells::Dict{Tuple{Int,Int},UInt8} = get_level_cells(z, level)
   
   # Trace the contour and return it
   sorted_contour::Vector{Vector{Vector{Float64}}} = trace_contour(x, y, z, level, level_cells, Vector{Float64})

   return sorted_contour
end

# The marching squares algorithm defines 16 cell types based on the edges that a 
# contour line enters and exits through. The edges of the cells are identified 
# using compass directions, while the vertices are ordered as follows:
#
#      N
#  4 +---+ 3
# W  |   |  E
#  1 +---+ 2
#      S
#
# Each cell type is identified with 4 bits, with each bit corresponding to a 
# vertex (MSB -> 4, LSB -> 1). A bit is set for vertex v_i is set if z(v_i) > h. 
# So a cell where a contour line only enters from the W edge and exits through 
# the N edge will have the cell type: 0b0111 Note that there are two cases where 
# there are two lines crossing through the same cell: 0b0101, 0b1010.
const N, S, E, W = (UInt8(1)), (UInt8(2)), (UInt8(4)), (UInt8(8))
const NS, NE, NW = N|S, N|E, N|W
const SN, SE, SW = S|N, S|E, S|W
const EN, ES, EW = E|N, E|S, E|W
const WN, WS, WE = W|N, W|S, W|E
const NWSE = NW | 0x10 # special (ambiguous case)
const NESW = NE | 0x10 # special (ambiguous case)

# Maps cell type to crossing types for non-ambiguous cells
const edge_LUT = (SW, SE, EW, NE, 0x0, NS, NW, NW, NS, 0x0, NE, EW, SE, SW)

# N, S, E, W
const next_map = ((0,1), (0,-1), (1,0), (-1,0))
const next_edge = (S,N,W,E)

@inline function _get_case(z, h)
   case = z[1] > h ? 0x01 : 0x00
   z[2] > h && (case |= 0x02)
   z[3] > h && (case |= 0x04)
   z[4] > h && (case |= 0x08)
   case
end

function get_level_cells(z, h::Number)
   cells = Dict{Tuple{Int,Int},UInt8}()
   x_ax, y_ax = axes(z)

   @inbounds for xi in first(x_ax):last(x_ax)-1
      for yi in first(y_ax):last(y_ax)-1
         elts = (z[xi, yi], z[xi + 1, yi], z[xi + 1, yi + 1], z[xi, yi + 1])
         case = _get_case(elts, h)

         # Contour does not go through these cells
         if iszero(case) || case == 0x0f
            continue
         end

         # Process ambiguous cells (case 5 and 10) using a bilinear interpolation
         # of the cell-center value.
         if case == 0x05
            cells[(xi, yi)] = 0.25*sum(elts) >= h ? NWSE : NESW
         elseif case == 0x0a
            cells[(xi, yi)] = 0.25*sum(elts) >= h ? NESW : NWSE
         else
            cells[(xi, yi)] = edge_LUT[case]
         end
      end
   end
   return cells
end

@inline function get_first_crossing(cell)
   if cell == NWSE
      return NW
   elseif cell == NESW
      return NE
   else
      return cell
   end
end

@inline function advance_edge(ind, edge)
   n = trailing_zeros(edge) + 1
   nt = ind .+ next_map[n]
   return nt, next_edge[n]
end

function get_next_edge!(cells::Dict, key, entry_edge::UInt8)
   cell = pop!(cells, key)
   if cell == NWSE
      if entry_edge == N || entry_edge == W
         cells[key] = SE
         cell = NW
      else #SE
         cells[key] = NW
         cell = SE
      end
   elseif cell == NESW
      if entry_edge == N || entry_edge == E
         cells[key] = SW
         cell = NE
      else #SW
         cells[key] = NE
         cell = SW
      end
   end
   return cell ⊻ entry_edge
end

# Given a cell and a starting edge, we follow the contour line until we either
# hit the boundary of the input data, or we form a closed contour.
function chase!(cells, curve, x, y, z, h, start, entry_edge, xi_range, yi_range, ::Type{VT}) where VT
   ind = start

   # When the contour loops back to the starting cell, it is possible
   # for it to not intersect with itself.  This happens if the starting
   # cell contains a saddle-point. So a loop is only closed if the
   # contour returns to the starting edge of the starting cell
   loopback_edge = entry_edge

   @inbounds while true
      exit_edge = get_next_edge!(cells, ind, entry_edge)

      push!(curve, interpolate(x, y, z, h, ind, exit_edge, VT))

      ind, entry_edge = advance_edge(ind, exit_edge)

      !((ind[1], ind[2], entry_edge) != (start[1], start[2], loopback_edge) &&
         ind[2] ∈ yi_range && ind[1] ∈ xi_range) && break
   end
   return ind
end


function trace_contour(x, y, z, h::Real, cells::Dict, VT)
   # Initialize the output Vector containing Vector of Tuples 
   contours::Vector{Vector{VT}} = []

   x_ax, y_ax = axes(z)
   xi_range = first(x_ax):last(x_ax)-1
   yi_range = first(y_ax):last(y_ax)-1

   @inbounds while length(cells) > 0
      contour_arr = VT[]

      # Pick initial box
      ind, cell = first(cells)

      # Pick a starting edge
      crossing = get_first_crossing(cell)
      starting_edge = 0x01 << trailing_zeros(crossing)

      # Add the contour entry location for cell (xi_0,yi_0)
      push!(contour_arr, interpolate(x, y, z, h, ind, starting_edge, VT))

      # Start trace in forward direction
      ind_end = chase!(cells, contour_arr, x, y, z, h, ind, starting_edge, xi_range, yi_range, VT)

      if ind == ind_end
         push!(contours, contour_arr)
         continue
      end

      # If loops reaches here then a reverse trace is needed
      ind, starting_edge = advance_edge(ind, starting_edge)

      if ind[2] ∈ yi_range && ind[1] ∈ xi_range
         # Start trace in reverse direction
         chase!(cells, reverse!(contour_arr), x, y, z, h, ind, starting_edge, xi_range, yi_range, VT)
      end
      push!(contours, contour_arr)
   end
   return contours
end

function interpolate(x::AbstractMatrix, y::AbstractMatrix, z::AbstractMatrix, h::Number, ind, edge::UInt8, VT)
   xi, yi = ind
   @inbounds if edge == W
      Δ = [y[xi,  yi+1] - y[xi,  yi  ], x[xi,  yi+1] - x[xi,  yi  ]].*(h - z[xi,  yi  ])/(z[xi,  yi+1] - z[xi,  yi  ])
      y_interp = y[xi,yi] + Δ[1]
      x_interp = x[xi,yi] + Δ[2]
   elseif edge == E
      Δ = [y[xi+1,yi+1] - y[xi+1,yi  ], x[xi+1,yi+1] - x[xi+1,yi  ]].*(h - z[xi+1,yi  ])/(z[xi+1,yi+1] - z[xi+1,yi  ])
      y_interp = y[xi+1,yi] + Δ[1]
      x_interp = x[xi+1,yi] + Δ[2]
   elseif edge == N
      Δ = [y[xi+1,yi+1] - y[xi,  yi+1], x[xi+1,yi+1] - x[xi,  yi+1]].*(h - z[xi,  yi+1])/(z[xi+1,yi+1] - z[xi,  yi+1])
      y_interp = y[xi,yi+1] + Δ[1]
      x_interp = x[xi,yi+1] + Δ[2]
   elseif edge == S
      Δ = [y[xi+1,yi  ] - y[xi,  yi  ], x[xi+1,yi  ] - x[xi,  yi  ]].*(h - z[xi,  yi  ])/(z[xi+1,yi  ] - z[xi,  yi  ])
      y_interp = y[xi,yi] + Δ[1]
      x_interp = x[xi,yi] + Δ[2]
   end
   return VT([x_interp, y_interp])
end

end