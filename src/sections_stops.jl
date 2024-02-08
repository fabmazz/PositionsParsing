function add_points_poly!(poly::AbstractDataFrame, points::AbstractDataFrame)
    for p in eachrow(points)
        ddi = haversine_d.(p.lat,p.lon,poly[!,:lat], poly[!,:lon])
        iis = sortperm(ddi)[1:2] #.+ (idx_p-1)
        i_att = minimum(iis)
        #println("$(ddi[iis])")
        if (ddi[iis[1]] < 5) || ((i_att == 1) && (2 in iis))
            ## non inserire il punto, assegna invece la fermata al punto iis[1]
            isel = iis[1]
            #println("Increment from $isel")
            poly[isel:end,:stopIdx] .+=1
        else
            ## insert the point after i_att
            r = poly[i_att,:]
            insert!(poly,i_att+1,merge(r,
                            (lat=p.lat,lon=p.lon,) ))
            #println("Insert point and increment at $(i_att+1)")
            ## increment
            poly[i_att+1:end,:stopIdx] .+= 1
        end
    end
end

function find_good_sections(poly::AbstractDataFrame, SEC_LEN_IDEAL::Real=500, lim_large::Real=200)
    istops = [1; findall(@. (poly.stopIdx[2:end]-poly.stopIdx[1:end-1])>0); length(poly)]
    #idx_step = zeros(Int,size(istops))
    ip = 1
    dists=zeros(length(istops)-1)
    for (ii,x) in enumerate(istops[2:end])
        dists[ii] = sum(distance_df(poly, ip, x))
        ip =x
    end
    ns=length(istops)
    sections = DataFrame(Dict("stop_s"=>collect(1:ns-1), "stop_e"=>collect(2:ns),"dist"=>dists));
    secs = sections #copy(sections);

    nn = ns-1
    for i=1:nn
        if i > size(secs,1)-1
            ##reached the end
            break
        end
        sdist = secs[i,:dist]
        if sdist> SEC_LEN_IDEAL
            continue
        elseif (sdist < SEC_LEN_IDEAL) && (secs[i+1,:dist]< SEC_LEN_IDEAL) &&(secs[i+1,:dist]+sdist < SEC_LEN_IDEAL+lim_large)
            ## fuse
            secs[i,:dist] += secs[i+1,:dist]
            secs[i,:stop_e] = secs[i+1,:stop_e]
            #println("Delete row $(i+1)")
            deleteat!(secs,i+1)
        end
    end
    secs
end