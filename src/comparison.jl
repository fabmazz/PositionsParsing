function unify_patterns_from_stops(patterns_sel::Dict{String,Any})
    nstops = 0
    key_sel = Set(keys(patterns_sel))
    for (k,dd) in patterns_sel
        #println("key $k")
        st1 = map(Stop,dd["stops"])
        if !(k in key_sel)
            continue
        end
        contained = false
        for k2 in key_sel
            if k2==k
                continue
            end
            p2=patterns_sel[k2]
            #println(p2)
            st2 = map(Stop,p2["stops"])
            mk=isnothing.(indexin(st1,st2))
            contains_all = !any(mk)
    
            if contains_all
                #println(" $k is contained in $k2, removing from good ones")
                contained = true
                break
            else
                l=length(st1[mk])
                #println("Stops not found in $k2: $(l)")
                if l==1
                    println(st1[mk])
                end
            end
        end
        if contained
            delete!(key_sel,k)
        end
    end
    key_sel
end

function get_longest_pattern_dist(polys::Dict{S,DataFrame}) where S<:AbstractString
    longest=""
    dist_long = 0.0

    for (k,pl) in polys
        #ax.plot(pl.lon,pl.lat,"o--", label=k)
        
        dd = sum(distance_df(pl))
        #println("distance tot: $(dd)")
        if dd > dist_long
            dist_long = dd
            longest=k
        end
        #plt.plot(pl.lon,pl.lat,"o--",alpha=0.5)
    end

    dlong = polys[longest]
    keep_poly = Dict(k => true for k in keys(polys))
    
    for (k,pl) in polys
        if k == longest
            continue
        end
        p1=pl
        print("$k ")
        #println(size(p1))
        ## calculate distance between points
        dd=haversine_d.(p1.lon', p1.lat', dlong.lon, dlong.lat)
        ii=argmin(dd,dims=1)
        #println(ii)
        idcs_2 = vec(map(x->x[1],ii))
        dist = vec(dd[ii])
        dist_m = (round.(Int,dist))

        idc_far=findall(dist_m .> 20)
        print(idc_far,"   ")
        accept = false
        if length(idc_far) < 1
            ## do nothing
        elseif length(idc_far) < 2
            il = idc_far[1]
            diff_dist = sum(distance_df(p1[il-1:il+1,:])) - sum( distance_df(dlong[idcs_2[il-1]:idcs_2[il+1],:]) )
            if diff_dist < 10
                println("REMOVE")
                accept = true
            else
                println("KEEP")
                accept=false
            end
        else
            ## we have two or more
            is=idc_far[1]
            ic = idc_far[1]
            #accept = true
            for (ik, idc) in enumerate(idc_far)
                if ik == 1
                    continue
                elseif ik > idc_far[ik-1]+1
                    ### break of streak
                    ## the idc is [is:ik-1]
                    println("Section between $is and $(ik-1)")
                    ie = ik-1
                    diff_dist = sum(distance_df(p1[is-1:ie+1,:])) - sum( distance_df(dlong[idcs_2[is-1]:idcs_2[ie+1],:]) )
                    if diff_dist > 10
                        println("KEEP")
                        accept=false
                    else
                        accept = true
                    end
                    is = ik
                end
            end
        end
        

        if accept
            keep_poly[k] = false
        end
    end

    [p[1] for p in filter(p -> p[2], keep_poly)]
end