#=function cut_trace_time(df_mat, mins):
    tdiff=df_nol["timerec"].diff().dt.seconds()
    tdiff[0]=0
    iidcs=np.where(tdiff > mins*60)[0]

    s=0
    datfs = []
    for i in iidcs:
        datfs.append(df_nol[s:i])
        s=i
    end
    datfs.append(df_nol[s:])
    return datfs
end
=#

function filter_positions_bbox(df, BBOX_LIMS)
    valid =@. (df.lat > BBOX_LIMS[1][1]) & (df.lat < BBOX_LIMS[1][2]) & (df.lon > BBOX_LIMS[2][1]) & (df.lon < BBOX_LIMS[2][2]);
    df_corr = df[valid,:]
    df_corr
end


function cut_trace_time(trace, mins_pass::Number)
    #distance = haversine_d.(trace.lat[1:end-1],trace.lon[1:end-1],trace.lat[2:end],trace.lon[2:end])
    nt = length(trace.timerec)
    times = trace.timerec[2:end].-trace.timerec[1:end-1]
    idcs = findall(times .> mins_pass*60)
    newidcs = idcs.+1
    if length(idcs) == 0
        return [trace]
    end

    out_ls = DataFrame[]
    si = 1
    for ind in idcs
        push!(out_ls,trace[si:ind,:])
        si = ind+1
    end
    finidc = idcs[end]
    if finidc < nt
        push!(out_ls, trace[finidc+1:end,:])
    end
    #println("Have $(length(out_ls)) dfs")
    out_ls
end

distance_df = tr -> haversine_d.(tr.lat[1:end-1],tr.lon[1:end-1],tr.lat[2:end],tr.lon[2:end])

#distance_df(tr::AbstractDataFrame) =  haversine_d.(tr.lat[1:end-1],tr.lon[1:end-1],tr.lat[2:end],tr.lon[2:end])
_distance_df(mdf::AbstractDataFrame,i0::Int,i1::Int) = haversine_d.(mdf.lat[i0:i1-1],mdf.lon[i0:i1-1],mdf.lat[i0+1:i1],mdf.lon[i0+1:i1])

function row_avg(df, j, i, count)
    for key in [:lat, :lon, :heading, :speed]
        if ((key==:heading) || (key==:speed)) && ismissing(df[j,key])
            df[j,key] = df[i,key] ## ignore None heading
        else
            df[j,key] = (df[j,key]*count[j] + df[i,key])/(count[j]+1)
        end
    end
    count[j]+=1
    count[i]=0
end

function average_pos_equaltime(spec_df)
    lendf = size(spec_df,1)
    count = ones(Int,lendf) #np.ones(len(spec_df),dtype=np.int_)
    secs_diff = spec_df.timerec[2:end].-spec_df.timerec[1:end-1];
    trec = spec_df[!,:timerec]

    idcs_eq = findall(secs_diff.==0) #np.where(secs_diff == 0)[0]
    if length(idcs_eq) == 0
        return spec_df
    end
    for i in idcs_eq
        #@assert trec[i] == trec[i+1]
        for j=i:-1:1
            if count[j] > 0
                row_avg(spec_df, j, i+1, count)
                break
            end
        end
    end
    #println("Have $(length(idcs_eq)) value with same time: $count")
    return spec_df[count.>0,:]
end

function _calc_new_speed(lat, lon, timerec, i_next,i_prev)
    dist_n = haversine_d(lat[i_prev],lon[i_prev], lat[i_next],lon[i_next])
    if dist_n == 0
        dist_n = 0.0001
    end
    tdiff = timerec[i_next]-timerec[i_prev]
    speed_new = dist_n*3.6/tdiff
    return dist_n, tdiff, speed_new
end

function remove_jumps_df(df::AbstractDataFrame, dist_all::Vector{T}, threshold::N = 100, debug::Bool=false) where T <: AbstractFloat where N <: Number
    if size(df,1) < 3
        ## with 2 points there are 1 arc with speed
        return df
    end

    timerec = df[!,:timerec] #.- df[1,:timerec]
    diff_sec = df[2:end,:timerec] .- df[1:(end-1),:timerec]
    lat = df[!,:lat]
    lon = df[!,:lon]
    if any(diff_sec .== 0)
        throw(ArgumentError("Some time differences are 0. Refusing to proceed"))
    end
    ##get speed per change
    speed = @. dist_all *3.6 / diff_sec
    #println("speeds: $speed")
    ## check the limit
    highspeed = (speed .>= threshold)
    if !any(highspeed)
        return df
    end
    idcs_hi = findall(highspeed)
    #println("is too high: $idcs_hi")
    keep = ones(Bool, size(df,1))
    for i in idcs_hi
        i_p = i
        i_2 = i+1
        if (i==1) 
            if (!highspeed[i+1]) 
                ## means that the first jump is far too fast, but not the following
                ## could be counteracted by the second jump
                keep[i] = false
            end
        elseif highspeed[i-1]
            ## previous jump was too fast also
            # calc new speed
            dist_n, tdiff, speed_new = _calc_new_speed(lat, lon, timerec, i+1, i-1)
            #println("from $(i-1) to $(i+1)")

            remove = ((dist_all[i]/dist_n) > 10) && (dist_n < 500) && !(speed_new > threshold)
            
            if remove
                keep[i] = false
                if debug
                    println(f"ratio {(dist_all[i]/dist_n):4.1f} , dist new {dist_n:4.1f} , speed new {speed_new:4.1f} remove {remove}")
                end
            end
        end
        ## new condition, after the previous jump
        if keep[i] && (i>2) && (highspeed[i-2]) 
            ## the previous jump was not too fast, but the second one was
            # calc new speed
            dist_n, tdiff, speed_new = _calc_new_speed(lat, lon, timerec, i+1, i-2)
            #println("from $(i-1) to $(i+1)")

            remove = ((dist_all[i]/dist_n) > 50) && (dist_n < 500) && !(speed_new > threshold)
            
            if remove
                keep[i] = false
                keep[i-1] = false
                #println("Remove $i and $(i-1)")
                if debug
                    println(f"ratio {(dist_all[i]/dist_n):4.1f} , dist new {dist_n:4.1f} , speed new {speed_new:4.1f}, remove 2-jump")
                end
            end
        end
    end
    df[keep,:]
end
_speed(dist, time) = @. dist *3.6 / time

function _check_cut_at_index!(mdf::AbstractDataFrame,i::Int, jend::Int, keep::Vector,speed_thrs::Number)
    #jend = jend-1
    dist = haversine_d.(mdf.lat[i],mdf.lon[i], mdf.lat[i+1:jend],mdf.lon[i+1:jend])
    ibest = argmin(dist)#+i
    #println("$i $jend, dist: $dist, min at $(argmin(dist)) -> $(ibest+i)")
    mindist = dist[ibest] == 0 ? 0.00001 : dist[ibest]
    jbest = ibest+i
    #disttot = sum((@. haversine_d(mdf.lat[i:jbest-1],mdf.lon[i:jbest-1],mdf.lat[i+1:jbest],mdf.lon[i+1:jbest]) ))#MatoParsing.distance_df(mdf[i:jbest,:]))
    disttot = sum(_distance_df(mdf, i, jbest))
    #println("distlong: $disttot, mindist: $mindist,  ratio: $(disttot/mindist)")
    tdiff = mdf.timerec[jbest]-mdf.timerec[i]
    ## new speed
    newspeed = _speed(mindist, tdiff) #mindist *3.6 / tdiff
    oldspeed = _speed(disttot, tdiff)
    
    if (disttot/mindist) > 15 && (newspeed < speed_thrs)
        ## cut trace
        keep[i+1:jbest-1] .= false
        println(f"Remove jump of {disttot:4.1f} m -> {mindist:3.1f}, old speed: {oldspeed:3.1f} new speed: {newspeed:3.1f}")
        if tdiff > 60
            println("Warning: removing loop with wrong speed creates time delta of $tdiff s")
        end
    end
end

function remove_large_jumps_df(mdf::AbstractDataFrame, dist::Vector{T}, speed_thrs::Number) where T<:AbstractFloat

    #dist = MatoParsing.distance_df(mdf)
    lendf = size(mdf,1)
    times = mdf.timerec[2:end] .- mdf.timerec[1:end-1]
    speed = _speed(dist,times) #@. dist *3.6 / times;

    mk = speed .> speed_thrs
    if !any(mk) || all(mk)
        ## in both cases there is nothing that can be done
        return mdf
    end
    iis=findall(mk)

    keep = ones(Bool, size(mdf,1))
    for i in iis
        #gs = cumsum(mk[i:end])
        if !keep[i]
            continue
        end
        #println(gs)
        jend=-2
        for j=i:length(mk)
            if !mk[j]
                #push!(jend,j)
                jend = j
                break
            end
        end
        if jend==-2 
            #println("Reached end without finding the element")
            #print(mk)
            ## remove the last part
            jcut = lendf
            disttot = sum(_distance_df(mdf, i, jcut))
            tdiff = mdf.timerec[jcut]-mdf.timerec[i]
            avg_speed = _speed(disttot, tdiff)
            if (avg_speed > speed_thrs) && (tdiff < 60)
                keep[i+1:jcut] .= false
                println(f"Remove sec with speed {avg_speed:3.1f} and dist {disttot:3.1f}")
            end
        else
            _check_cut_at_index!(mdf, i, jend, keep, speed_thrs)
        end
        
    end

    mdf[keep,:]
end