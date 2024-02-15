using Polyline

using LibGEOS

function get_pattern_polyline(patterns::AbstractDict, code::String)
    pl=patterns[code]["patternGeometry"]["points"]

    tracePoly = DataFrame(decodePolyline(pl),["lat","lon"])
    tracePoly[!,:idx] = collect(1:size(tracePoly,1))

    tracePoly
end

function get_pattern_stops(patterns::AbstractDict, code::String)
    st = patterns[code]["stops"]
    df = DataFrame(st)
    df[!,:stopIdx]=collect(1:size(df,1))
    
    return df
end

function add_midpoints_idc!(geom::AbstractDataFrame, npoints::Integer, idc_point::Integer; verbose=false)

    steplon = (geom[idc_point + 1, :lon] - geom[idc_point, :lon]) / (npoints + 1)
    steplat = (geom[idc_point + 1, :lat] - geom[idc_point, :lat]) / (npoints + 1)


    if verbose
        println("steps: ", (steplat,steplon))
    end
    
    ## adjust indices
    geom[idc_point+1:end,:idx] .+= npoints
    r = geom[idc_point,:]
    for i in 1:npoints
        newrow = (lat=geom[idc_point, :lat] +steplat * (i),
                lon=geom[idc_point, :lon] + steplon * (i), idx=idc_point+i)

        if verbose
            println(newrow)
        end

        #insert!(geom2, idc_point+i, [newp], names(geom))
        insert!(geom,idc_point+i, merge(r,newrow))
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

function fill_polyline_points(tracePoly::AbstractDataFrame, distMet::Number=20)
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
        dist_poly = distance_df(polyLong)
    end
    @assert sum(dist_poly.==0) == 0
    polyLong[!,:idx] = collect(1:size(polyLong,1))
    #@assert sum(dist_poly.>20) == 0
    polyLong
end

function set_sections_poly!(polyLong::AbstractDataFrame, sec_length::Number = 500,dist_start_end::Number = 30)

    ## assign sections
    dist_poly = distance_df(polyLong)
    dist_cum = cumsum([0.0;dist_poly])
    i_s = findfirst(dist_cum.> dist_start_end)
    ilast=findlast(@. (dist_cum[end] - dist_cum) > dist_start_end)
    print("idcs: $i_s, $ilast")
    secIdxAll = ones(Int,length(dist_cum))
    secIdxAll[1:i_s-1] .= 1
    iep= ilast-1
    secIdxAll[i_s:iep] = (floor.(cumsum(dist_poly[i_s:iep])./ sec_length)).+2
    secIdxAll[ilast:end] .= secIdxAll[iep]+1

    polyLong[!,:section] = secIdxAll;
    
end

function calc_dist_sections(polyLong::AbstractDataFrame)
    idc_sec_change=findall(Bool.(polyLong.section[2:end]-polyLong.section[1:end-1])).+1;

    #idc_sel2 ::Vector{Int64}= [1;idc_sec_change;size(polyLong,1)]
    #idc_sec_change

    secsPoly=polyLong[idc_sec_change,:]

    r=DataFrame( (dist=distance_df(polyLong),section=polyLong[1:end-1,:section]))
    dist_sec = combine(groupby(r,:section),:dist => sum)
    dist_sec
end

### closest point of the polyline to the trace

function find_closest_poly_trace(trace::AbstractDataFrame, poly::AbstractDataFrame, dist_accept::Number=50)
    dist_mat =haversine_d.(trace[!,:lat],trace[!,:lon], poly[!,:lat]', poly[!,:lon]');

    cartidx_nearest=argmin(dist_mat, dims=2)
    dist_actual=dist_mat[cartidx_nearest];
    accept = vec(dist_actual .< dist_accept)

    idx_nearest = vec(map(x->x[2],cartidx_nearest))
    #trace[!,:polyPoint] = idx_nearest 
    #trace[!,:secIdx]=
    sec_nearest = poly[idx_nearest,:section];
    return accept, idx_nearest, sec_nearest, dist_mat
end


function find_closest_poly_trace!(trace::AbstractDataFrame, poly::AbstractDataFrame)
    idx_nearest, sec_nearest,_ = find_closest_poly_trace(trace, poly)
    trace[!,:polyPoint] = idx_nearest
    trace[!,:secIdx] = sec_nearest

end


function remap_df_points(points)
    lats=Float64[]
    lons=Float64[]
    for p in points
        push!(lats,LibGEOS.getGeomX(p))
        push!(lons,LibGEOS.getGeomY(p))
    end
    DataFrame(Dict(:lat=>lats,:lon=>lons))
end

function get_nearest_point_vecs(trace::AbstractDataFrame, polySections::AbstractDataFrame; addsections=true)
    trace_libgeos = LibGEOS.LineString(map(x -> [x.lat,x.lon], eachrow(trace)))
    points_closest = map(x->LibGEOS.nearestPoints(trace_libgeos,LibGEOS.Point(x.lat,x.lon))[1], eachrow(polySections));
    newpoints= remap_df_points(points_closest)
    if addsections
        newpoints[!,:secIdx] = polySections.section
    end
    newpoints
end



function fix_reverting_section_postrace!(trace::AbstractDataFrame, polyLong::AbstractDataFrame,dist_thresh::Number=10)

    pos_sect = trace.secIdx
    len_trace=size(trace,1)
    diff = pos_sect[2:end] - pos_sect[1:end-1]
    ii = findall(diff .< 0)
    if length(ii) == 0
        return
    end
    for i in ii
        #change from point from point i to point i+1
        println("i $i")
        if diff[i]<-1
            ## looping around (like with the tram 9 at the end) 
            trace[i+1,:secIdx] = trace[i,:secIdx]
            ## fix polyPoint to closest to the section 
            secsPoints = polyLong[polyLong.section.== trace[i,:secIdx],:]
            dist = haversine_d.(secsPoints.lat, secsPoints.lon, trace[i+1,:lat], trace[i+1,:lon])
            closest = secsPoints[argmin(dist),:]
            trace[i+1, :polyPoint] = closest.idx[1]

            ## recalculate diff
            diff = trace.secIdx[2:end] - trace.secIdx[1:end-1]
            println("Fixed point going too back in the section idx")
        elseif (i==len_trace-1) && (i>2)
            ## reverting in the last idx
            ifix = i+1 ## last 
            @assert ifix == len_trace
            if (haversine_d(trace, i, i-1) < dist_thresh) && (haversine_d(trace,i,i+1) < dist_thresh)
                ## 
                trace[i,:secIdx] =  trace[i+1,:secIdx]
            elseif (haversine_d(trace,i,i+1) < dist_thresh)
                trace[i+1,:secIdx] = trace[i,:secIdx]
            end

        elseif (i > 1) && diff[i-1] == 1 && diff[i+1]==1
            ## i has the section +1
            # fix section of prev point
            ifix0 = i+1
            ifix1 = i-1
            dprev = haversine_d(trace, i-1, i)
            dnext =  haversine_d(trace, i,i+1)
            polyPoints = trace.polyPoint

            if (dprev < dist_thresh) && (dprev < dnext)
                ## the point i should have the previous section, eg: [17, 18 (i), 17, 18] -> [17, 17(i), 17, 18]
                trace[i,:secIdx] = trace[i-1,:secIdx]
                ## set the polyPoint to educated guess
                trace[i,:polyPoint] = polyPoints[i-1] + Int(round(Int, (polyPoints[i+1]-polyPoints[i-1])/2))
                println("Changed secIdx at $ifix0 ")
            elseif (dnext < dist_thresh) && (dnext < dprev)
                ## the point i+1 should have the previous section, eg: [17, 18 (i), 17, 18] -> [17, 18(i), 18, 18]
                trace[i+1,:secIdx] = trace[i,:secIdx]
                ## set the polyPoint to educated guess
                trace[i+1,:polyPoint] = polyPoints[i] + Int(round(Int, (polyPoints[i+2]-polyPoints[i])/2))
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
        if !(length(idcs_df)==1)
            #println("Find single sect, length dist of section $sidx is $(length(idcs_df)), mask has $(sum(mask)) values. This should not happen")
            #println(res)

            continue
        end
        
        mSec = secsPoly[secsPoly.section.==sidx+1,:]
        if length(mSec) > 0
            println("secfix sec: $sidx, $(sum(mask)); $(size(mSec)), $idcs_df")
            dist =  haversine_d.(trace[mask,:lat],trace[mask,:lon],mSec.lat, mSec.lon)
            
            @assert length(dist)==1

            ii=idcs_df[1]
            if dist[1] < dist_thresh
                trace[ii,:secIdx] = sidx+1
                changed[1] = true
                ## fix polyPoint
                trace[ii,:polyPoint] = mSec.idx[1]
            end
            _= any(changed) ? println("Raised secIdx at $ii, previous: $sidx ") : nothing
        end

        if sidx > 1 && !changed[1]
            ## do the same with index sidx -1
            mSec = secsPoly[secsPoly.section.==sidx-1,:]
            if length(mSec) > 0
                println("sec: $sidx, $(sum(mask)); $(size(mSec)), $idcs_df")
                dist =  haversine_d.(trace[mask,:lat],trace[mask,:lon],mSec.lat, mSec.lon)
                @assert length(dist)==1
    
                ii=idcs_df[1]
                if dist[1] < dist_thresh
                    trace[ii,:secIdx] = sidx-1
                    changed[ii] = true
                    ## fix polyPoint
                    trace[ii,:polyPoint] = mSec.idx[1]
                end
                _= any(changed) ? println("Raised secIdx at $ii, previous: $sidx ") : nothing
            end
        end

    end
    
end
        

function add_sec_points_postrace(trace::AbstractDataFrame, secsPoly::AbstractDataFrame, newpoints, dist_thresh::Number=20; 
    timecol::Symbol=:timerec, verbose=false)
    df = copy(trace)
    nump = size(newpoints,1)

    for (i,p) in enumerate(eachrow(newpoints))
        if p.secIdx == 1
            continue
        end
        g = verbose ? println("Point $p") : nothing
        #p = newpoints[i,:]
        isec = findfirst(secsPoly.section.==p.secIdx)
        idc_pol = secsPoly[isec,:idx]
        sect_idx = secsPoly[isec,:section]
        
        idf = findfirst(df.polyPoint .>=idc_pol)
        if isnothing(idf)
            ## no point with higher
            #println("No point")
            continue
        elseif df.secIdx[idf] > sect_idx#isnothing(idf)
            ## index is not present
            #println("section idx is too large, sec $sect_idx not present?")
            continue
        end
        if idf == 1
            ## first point in a section
            println("First point")
            continue
        end
        #ip = i_n-1
        #dist = MatoParsing._distance_df(df,idf-1,idf)
        g=verbose ? println("i $i -> sec $(sect_idx) idf $idf idc_pol: $idc_pol") : nothing
        if df.polyPoint[idf] == idc_pol
            ## close to same point
            dd  =haversine_d(df.lat[idf],df.lon[idf],secsPoly.lat[isec], secsPoly.lon[isec])
            if dd < dist_thresh
                ## do not add point
                println("Skip adding point sec $sect_idx, too close")
                continue
            end
        end
        tp = df[idf-1,timecol]
        tn = df[idf,timecol]
        
        if df.secIdx[idf] == df.secIdx[idf-1]
            throw(AssertionError("Section indices at $idf, $(idf-1) are not compatible with the closest polypoints for section $(sect_idx):\n $(df[idf-1:idf,:])\n Are you sure you haven't already run this method?"))
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
            println("new point == previous for sec $(sect_idx), avoid inserting, idcPoly: prev=$(df.polyPoint[idf-1]), next=$(df.polyPoint[idf])")
            ## do not insert point, modify previous point
            df[idf-1,:secIdx] = sect_idx
            ins=false
        end
        if ins
            if verbose
                println(f"Insert new point at {p.lat:5.4f} {p.lon:5.4f} d0 {dist_prev:3.1f} d1 {dist_next:3.1f} t_new {tnew_i} idcPoly {idc_pol}")
            end
            row=(df[idf-1,:])
            ## insert row in the DataFrame
            add =(; lat=p.lat,lon=p.lon, timecol => tnew_i, polyPoint=idc_pol,secIdx=sect_idx) 
            newrow =merge(row, add)
            #timecol=tnew_i
            #println(add)
            insert!(df,idf,newrow)
        end

        
    end

    df
end

