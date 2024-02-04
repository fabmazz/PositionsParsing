function divide_trace_segment!(geom::AbstractDataFrame, dist, max_dist, idc_point; verbose=false)
    @assert size(geom,1) == length(dist) + 1

    npoints = Int(floor(dist[idc_point] / max_dist))

    steplon = (geom[idc_point + 1, :lon] - geom[idc_point, :lon]) / (npoints + 1)
    steplat = (geom[idc_point + 1, :lat] - geom[idc_point, :lat]) / (npoints + 1)


    if verbose
        println("steps: ", (steplat,steplon))
    end
    
    ## adjust indices
    geom[idc_point+1:end,:idx] .+= npoints
    for i in 1:npoints
        newrow = (lat=geom[idc_point, :lat] +steplat * (i + 1),
                lon=geom[idc_point, :lon] + steplon * (i + 1), idx=idc_point+i)

        if verbose
            println(newrow)
        end

        #insert!(geom2, idc_point+i, [newp], names(geom))
        insert!(geom,idc_point+i, newrow)
    end

    return geom, npoints
end

function add_midpoints_idc!(geom::AbstractDataFrame, npoints::Integer, idc_point::Integer; verbose=false)

    steplon = (geom[idc_point + 1, :lon] - geom[idc_point, :lon]) / (npoints + 1)
    steplat = (geom[idc_point + 1, :lat] - geom[idc_point, :lat]) / (npoints + 1)


    if verbose
        println("steps: ", (steplat,steplon))
    end
    
    ## adjust indices
    geom[idc_point+1:end,:idx] .+= npoints
    for i in 1:npoints
        newrow = (lat=geom[idc_point, :lat] +steplat * (i + 1),
                lon=geom[idc_point, :lon] + steplon * (i + 1), idx=idc_point+i)

        if verbose
            println(newrow)
        end

        #insert!(geom2, idc_point+i, [newp], names(geom))
        insert!(geom,idc_point+i, newrow)
    end
    
end

function add_midpoints_trace_all!(geom::AbstractDataFrame, max_dist::Number; verbose=false)
    
    dist = distance_df(geom)

    #geom2 = copy(geom)
    idc_added = 0
    npoints = @. dist / max_dist
    for i in findall(dist .> max_dist)
        
        idc_point = i+idc_added
        npp = Int(floor(npoints[i]))
        @assert npp > 0
        add_midpoints_idc!(geom, npp, idc_point; verbose=verbose)

        idc_added += npp
    end

    return geom, idc_added
end

function fill_polyline_points(tracePoly::AbstractDataFrame,distMet::Number=20, sec_length::Number = 500)
    polyLong = copy(tracePoly);
    add_midpoints_trace_all!(polyLong,distMet; verbose=false)
    dist_poly = distance_df(polyLong)
    ## remove points with distance 0
    mask = ones(Bool, size(polyLong,1))
    mask[findall(dist_poly.==0)].=false
    #println(mean(mask)," ",sum(mask))
    polyLong = polyLong[mask,:]
    dist_poly = dist_poly[dist_poly.>0]

    ## assign sections
    
    idx_sec = Int.(floor.(cumsum([0.0;dist_poly]) / sec_length))

    polyLong[!,:section] = idx_sec;
    polyLong[!,:idx] = collect(1:size(polyLong,1))
    polyLong
end

### closest point of the polyline to the trace

function closest_poly_trace(trace::AbstractDataFrame, poly::AbstractDataFrame)
    
end