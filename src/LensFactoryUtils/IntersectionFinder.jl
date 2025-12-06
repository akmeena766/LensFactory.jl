module IntersectionFinder

# Based on: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
# Python implementation: https://github.com/sukhbinder/intersection/tree/master/intersect

# Functions to export
export get_intersection

function get_intersection(x1::AbstractVector{<:Real}, y1::AbstractVector{<:Real}, 
                          x2::AbstractVector{<:Real}, y2::AbstractVector{<:Real})::Vector{NTuple{2,Real}}
    # Check if curve-1 is consistent
    if !(length(x1) > 1) || !(length(y1) > 1) || length(x1) == length(y1) || 
        throw(ArgumentError("Incompatible input axes for input vectors."))
    end
    
    # Check if curve-2 is consistent
    if !(length(x2) > 1) || !(length(y2) > 1) || length(x2) == length(y2) || 
        throw(ArgumentError("Incompatible input axes for input vectors."))
    end    

    # Number of line segments in both curves
    n1 = length(x1) - 1
    n2 = length(x2) - 1

    xy1 = hcat(x1, y1)
    xy2 = hcat(x2, y2)

    dxy1 = diff(xy1, dims = 1)
    dxy2 = diff(xy2, dims = 1)

    ijc = Array{Any}(undef,n2)
    min_x1 = mvmin(x1)
	max_x1 = mvmax(x1)
	min_y1 = mvmin(y1)
	max_y1 = mvmax(y1)
    for k in 1:n2
        k1 = k + 1
        ijc[k] = findall( (min_x1 .≤ max(x2[k], x2[k1])) .& 
                          (max_x1 .≥ min(x2[k], x2[k1])) .&
                          (min_y1 .≤ max(y2[k], y2[k1])) .& 
                          (max_y1 .≥ min(y2[k], y2[k1])) )
        # Second column
        ijc[k] = [[ij[1], k] for ij in ijc[k]]
    end
    ij = vcat(ijc...)
    i = [vec[1] for vec in ij]
    j = [vec[2] for vec in ij]
    
    n = length(i)
    T = zeros(4, n)
    AA = zeros(4, 4, n)
    AA[[1, 2], 3, :] .= -1
    AA[[3, 4], 4, :] .= -1
    AA[[1, 3], 1, :] .= dxy1[i, :]'
    AA[[2, 4], 2, :] .= dxy2[j, :]'
    B = -hcat(x1[i], x2[j], y1[i], y2[j])'

    for k in 1:n        
        try
            T[:, k] .= AA[:, :, k] \ B[:, k]
        catch err        
            T[:, k] .= Inf        
        end
    end
    in_range = ( (T[1, :] .≥ 0.0) .& (T[2, :] .≥ 0.0) .& (T[1, :] .≤ 1.0) .& (T[2, :] .≤ 1.0) )'
    
    xy0=[]
    for ll in axes(in_range,2)
        if in_range[ll]
            push!(xy0,(T[3, ll], T[4, ll]))
        end
    end
    
    # Find unique solutions
   xy0_unique::Vector{NTuple{2,Real}}=[]
   l::Int = 0
   for i in axes(xy0,1)
      l = 0
      for j in i+1:length(xy0)
         norm = sqrt( (xy0[i][1]-xy0[j][1])^2 + (xy0[i][2]-xy0[j][2])^2)
         if norm < 1E-12
               l = l + 1
         end
      end
      if l == 0
         push!(xy0_unique, xy0[i])
      end
   end
   return xy0_unique
end

@inline function mvmin(x)
    return minimum.( eachrow( hcat(x[1:end-1, :], x[2:end, :]) ) )
end

@inline function mvmax(x)
    return maximum.( eachrow( hcat(x[1:end-1, :], x[2:end, :]) ) )
end

end