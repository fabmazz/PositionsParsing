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
        newrow = (lat=geom[idc_point, :lat] +steplat * (i),
                lon=geom[idc_point, :lon] + steplon * (i), idx=idc_point+i)

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
    if verbose
        println("npoints: $(round.(npoints,digits=2))")
    end
    for i in findall(dist .> max_dist)

        idc_point = i+idc_added
        @assert haversine_d(geom,idc_point, idc_point+1) == dist[i]
        npp = Int(floor(npoints[i]))
        #prlat = geom.lat[idc_point+1]
        #prlon = geom.lon[idc_point+1]
        #prev_dist = haversine_d(geom, idc_point, idc_point+1)

        @assert npp > 0
        if verbose
            println("i $i ,distance: $(round(dist[i],digits=2)) ,n_points $npp")
        end
        # fix double points
        if round(Int,dist[i])==npp*round(Int,max_dist)
            npp -=1
            #println("Removed one point")
            if npp == 0
                continue
            end
        end

        add_midpoints_idc!(geom, npp, idc_point; verbose=false)
        #@assert (prlat == geom.lat[idc_point+npp+1]) && (prlon == geom.lon[idc_point+npp+1])
        ## debugging
        #=
        newdist = distance_df(geom,idc_point, idc_point+npp+1)
        if verbose && any(newdist.>= max_dist)
            println("Prev dist: $prev_dist, Wrong distances, new dist: $newdist")
        end
        =#
        idc_added += npp
    end

    return geom, idc_added
end

function fill_polyline_points(tracePoly::AbstractDataFrame, distMet::Number=20, sec_length::Number = 500)
    polyLong = copy(tracePoly);
    add_midpoints_trace_all!(polyLong,distMet; verbose=false)
    dist_poly = distance_df(polyLong)
    ## remove points with distance 0
    mask = ones(Bool, size(polyLong,1))
    ii=(findall(dist_poly.==0))
    if length(ii)>0
        mask[findall(dist_poly.==0)].=false
        #println(mean(mask)," ",sum(mask))
        polyLong = polyLong[mask,:]
        dist_poly = distance_df(polyLong) #dist_poly[dist_poly.>0]
    end
    @assert sum(dist_poly.==0) == 0
    #@assert sum(dist_poly.>20) == 0

    ## assign sections
    
    idx_sec = Int.(floor.(cumsum([0.0;dist_poly]) / sec_length)) .+1

    polyLong[!,:section] = idx_sec;
    polyLong[!,:idx] = collect(1:size(polyLong,1))
    polyLong
end

### closest point of the polyline to the trace

function find_closest_poly_trace(trace::AbstractDataFrame, poly::AbstractDataFrame)
    dist_mat =haversine_d.(trace[!,:lat],trace[!,:lon], poly[!,:lat]', poly[!,:lon]');

    gg=argmin(dist_mat, dims=2)
    dist_acc=dist_mat[gg];

    idx_nearest = vec(map(x->x[2],argmin(dist_mat, dims=2)))
    #trace[!,:polyPoint] = idx_nearest 
    #trace[!,:secIdx]=
    sec_nearest = poly[idx_nearest,:section];
    return idx_nearest, sec_nearest, dist_mat
end


function find_closest_poly_trace!(trace::AbstractDataFrame, poly::AbstractDataFrame)
    idx_nearest, sec_nearest,_ = find_closest_poly_trace(trace, poly)
    trace[!,:polyPoint] = idx_nearest
    trace[!,:secIdx] = sec_nearest

end

using LibGEOS

function remap_df_points(points)
    lats=Float64[]
    lons=Float64[]
    for p in points
        push!(lats,LibGEOS.getGeomX(p))
        push!(lons,LibGEOS.getGeomY(p))
    end
    DataFrame(Dict(:lat=>lats,:lon=>lons))
end

function get_nearest_point_vecs(trace::AbstractDataFrame, polySections::AbstractDataFrame)
    trace_libgeos = LibGEOS.LineString(map(x -> [x.lat,x.lon], eachrow(trace)))
    points_closest = map(x->LibGEOS.nearestPoints(trace_libgeos,LibGEOS.Point(x.lat,x.lon))[1], eachrow(polySections));
    newpoints= remap_df_points(points_closest)
end



function fix_reverting_section_postrace!(trace::AbstractDataFrame,dist_thresh::Number=20)

    pos_sect = trace.secIdx
    diff = pos_sect[2:end] - pos_sect[1:end-1]
    ii = findall(diff .< 0)
    if length(ii) == 0
        return
    end
    for i in ii
        #seci = trace.sexIdx[i] ## previous
        println("i $i")
        if diff[i-1] == 1 && diff[i+1]==1
            ## i has the section +1
            # fix section of prev point
            ifix0 = i+1
            ifix1 = i-1
            if (haversine_d(trace, ifix0, ifix0-1) < dist_thresh)
                ## fix section
                trace[ifix0,:secIdx] = trace[ifix0+1,:secIdx]
                println("Changed secIdx at $ifix0 ")
            elseif (ifix1 > 1) &&  (haversine_d(trace, ifix1+1, ifix1) < dist_thresh)
                trace[ifix1,:secIdx] = trace[ifix1+1,:secIdx]
                println("Changed secIdx at $ifix1 ")
            end
        end
    end
end
function fix_section_pos_trace!(trace::AbstractDataFrame, secsPoly::AbstractDataFrame, dist_thresh::Number=11)

    ## find all sections with only 1 element
    res= combine( groupby(trace, :secIdx), :secIdx => length)
    secIdx1 = res[res.secIdx_length.==1,:secIdx]
    for sidx in secIdx1
        mask=(trace.secIdx.==sidx)
        idcs_df = findall(mask)  #df_fi = trace[mask,:]
        changed = zeros(Bool, length(idcs_df))
        mSec = secsPoly[secsPoly.section.==sidx+1,:]
        println(size(mSec), idcs_df)
        dist =  haversine_d.(trace[mask,:lat],trace[mask,:lon],mSec.lat, mSec.lon)
        print(length(dist))
        #iact =  dist_ .< dist_thresh
        #df_fi[iact,:secIdx] .= sidx+1
        for i in idcs_df
            if dist[i] < dist_thresh
                trace[i,:secIdx] = sidx+1
                changed[i] = true
            end
        end
        println("Raised secIdx in $(sum(changed)) cases")

        ##TODO if needed: assign the previous section index

    end
    
end
        

function add_sec_points_postrace(trace::AbstractDataFrame, secsPoly::AbstractDataFrame, newpoints, dist_thresh::Number=20)
    df = copy(trace)

    for i=1:size(newpoints,1)
        p = newpoints[i,:]
        idc_pol = secsPoly[i,:idx]
        sect_idx = secsPoly[i,:section]
        
        idf = findfirst(df.polyPoint .>=idc_pol)
        if isnothing(idf)
            ## no point with higher
            println("No point")
            continue
        elseif df.secIdx[idf] > sect_idx#isnothing(idf)
            ## index is not present
            println("section idx is too large, sec $sect_idx not present?")
            continue
        end
        if idf == 1
            ## first point in a section
            println("First point")
            continue
        end
        #ip = i_n-1
        #dist = MatoParsing._distance_df(df,idf-1,idf)
        println("i $i -> sec $(i+1) idf $idf idc_pol: $idc_pol")
        if df.polyPoint[idf] == idc_pol
            ## close to same point
            dd  =haversine_d(df.lat[idf],df.lon[idf],secsPoly.lat[i], secsPoly.lon[i])
            if dd < dist_thresh
                ## do not add point
                println("Skip adding point, too close")
                continue
            end
        end
        tp = df.timerec[idf-1]
        tn = df.timerec[idf]
        
        if df.secIdx[idf] == df.secIdx[idf-1]
            throw(AssertionError("Section indices are not compatible with the closest polypoints for section $(i+1). Are you sure you haven't already run this method?"))
        end
        dist_prev = haversine_d(df.lat[idf-1],df.lon[idf-1], p.lat, p.lon)
        dist_next = haversine_d(p.lat,p.lon,df.lat[idf],df.lon[idf])
        tnew_i = round(Int,
            (dist_prev*tn+dist_next*tp)/(dist_prev+dist_next)
            )
        ins =true
        if (dist_next==0)
            println("Avoid inserting point of same section and 0 dist to next, sec $(i+1)")
            ins =false
        elseif (dist_prev == 0)
            println("new point == previous for sec $(i+1), avoid inserting, idcPoly: prev=$(df.polyPoint[idf-1]), next=$(df.polyPoint[idf])")
            ## do not insert point, modify previous point
            df[idf-1,:secIdx] = i+1
            ins=false
        end
        if ins
            println(f"Insert new point at {p.lat:5.4f} {p.lon:5.4f} d0 {dist_prev:3.1f} d1 {dist_next:3.1f} t_new {tnew_i} idcPoly {idc_pol}")
            row=(df[idf-1,:])
            ## insert row in the DataFrame
            insert!(df,idf,merge(row,
                        (lat=p.lat,lon=p.lon,timerec=tnew_i,polyPoint=idc_pol,secIdx=i+1) ))
        end

        
    end

    df
end