module AreaFinder

export shoelace

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

end