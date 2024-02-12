

function process_polyline_constant_sections(tracePoly::AbstractDataFrame, POLY_DIST::Integer=20, SEC_LENGTH::Integer=500,FIRST_LAST_SEC::Integer=30 )
    polyLong = fill_polyline_points(tracePoly,POLY_DIST)
    set_sections_poly!(polyLong, SEC_LENGTH, FIRST_LAST_SEC);

    ## find sections first element

    polyLong
end

function find_section_points_polyLong(polyLong::AbstractDataFrame)
    idc_sec_change=findall(Bool.(polyLong.section[2:end]-polyLong.section[1:end-1])).+1;
    idc_sec_change
    secsPoly=vcat(DataFrame(polyLong[1,:]),polyLong[idc_sec_change,:])
    secsPoly
end

"""
Process the polyline creating sections corresponding to the stops (stops need to be ordered)
"""
function process_polyline_sections_stops(poly::AbstractDataFrame, stops::AbstractDataFrame, POLY_DIST::Real=20, SEC_LEN_IDEAL::Real=500, DIST_EDGE_SEC::Real=30)
    points = get_nearest_point_vecs(poly, stops; addsections=false);

    pc =copy(poly);
    pc[!,:stopIdx] .= 0
    add_points_poly!(pc, points)

    secs=find_good_sections(pc, SEC_LEN_IDEAL)

    pc[:,:section] = map( si -> findfirst(@. (si >= secs.stop_s) & (si < secs.stop_e)),pc.stopIdx);
    
    ## add points
    
    polyLong = fill_polyline_points(pc,POLY_DIST)
    send = polyLong[end,:section]
    if size(polyLong[polyLong.section.==send,:],1) ==1
        println("Last section has only one element, removing")
        polyLong[end,:section] = polyLong[end-1,:section]
    end


    
    mdist = distance_df(polyLong)
    dist_cum = cumsum([0.0;mdist])
    i_s = findfirst(dist_cum.> DIST_EDGE_SEC)
    ilast=findlast(@. (dist_cum[end] - dist_cum) > DIST_EDGE_SEC)

    polyLong[i_s:end,:section].+=1
    polyLong[ilast:end,:section] .+=1

    polyLong
end