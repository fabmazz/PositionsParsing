
distance_df(tr::AbstractDataFrame) =  haversine_d.(tr.lat[1:end-1],tr.lon[1:end-1],tr.lat[2:end],tr.lon[2:end])
distance_df(mdf::AbstractDataFrame,i0::Int,i1::Int) = haversine_d.(mdf.lat[i0:i1-1],mdf.lon[i0:i1-1],mdf.lat[i0+1:i1],mdf.lon[i0+1:i1])

haversine_d(df::AbstractDataFrame, i::Integer, j::Integer) = haversine_d(df.lat[i],df.lon[i],df.lat[j],df.lon[j])

timediff(mdf::AbstractDataFrame) = mdf.timerec[2:end] .- mdf.timerec[1:end-1]
timediff(mdf::AbstractDataFrame,i0::Integer,i1::Integer) = mdf.timerec[i0+1:i1] .- mdf.timerec[i0:i1-1]

speeds_df(mdf::AbstractDataFrame) = _speed(distance_df(mdf), timediff(mdf))

function filter_positions_bbox(df, BBOX_LIMS)
    valid =@. (df.lat > BBOX_LIMS[1][1]) & (df.lat < BBOX_LIMS[1][2]) & (df.lon > BBOX_LIMS[2][1]) & (df.lon < BBOX_LIMS[2][2]);
    df_corr = df[valid,:]
    df_corr
end


function cut_trace_time(trace, mins_pass::Number, time_col::Symbol=:timerec)
    #distance = haversine_d.(trace.lat[1:end-1],trace.lon[1:end-1],trace.lat[2:end],trace.lon[2:end])
    nt = size(trace,1)
    times = trace[2:end,time_col].-trace[1:end-1,time_col]
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

function cut_trace_time_dist(trace::AbstractDataFrame, mins_pass::Number, dist_meters::Number)
    #distance = haversine_d.(trace.lat[1:end-1],trace.lon[1:end-1],trace.lat[2:end],trace.lon[2:end])
    nt = length(trace.timerec)
    times = trace.timerec[2:end].-trace.timerec[1:end-1]
    dist = distance_df(trace)
    idcs = findall(@. (times > mins_pass*60) | (dist > dist_meters))
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


function row_avg(df, j, i, count)
    for key in [:lat, :lon, :heading]# :speed]
        if ((key==:heading) || (key==:speed)) && ismissing(df[j,key])
            df[j,key] = df[i,key] ## ignore None heading
        else
            df[j,key] = (df[j,key]*count[j] + df[i,key])/(count[j]+1)
        end
    end
    count[j]+=1
    count[i]=0
end

function average_pos_equaltime(spec_df::AbstractDataFrame, time_col::Symbol=:timerec)
    lendf = size(spec_df,1)
    count = ones(Int,lendf) #np.ones(len(spec_df),dtype=np.int_)
    secs_diff = spec_df[2:end,time_col].-spec_df[1:end-1,time_col];
    #trec = spec_df[!,:timerec]

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

function remove_jumps_df(df::AbstractDataFrame, dist_all::Vector{T}, threshold::N = 100, debug::Bool=false; time_col::Symbol=:timerec) where T <: AbstractFloat where N <: Number
    if size(df,1) < 3
        ## with 2 points there are 1 arc with speed
        return df
    end

    timerec = df[!,time_col] #.- df[1,:timerec]
    diff_sec = df[2:end,time_col] .- df[1:(end-1),time_col]
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

function _check_cut_at_index!(mdf::AbstractDataFrame,i::Int, jend::Int, keep::Vector,speed_thrs::Number, time_col::Symbol)
    #jend = jend-1

    ## i is the index that is the last correct point
    dist = haversine_d.(mdf.lat[i],mdf.lon[i], mdf.lat[i+1:jend],mdf.lon[i+1:jend])
    ibest = argmin(dist)#+i
    #println("$i $jend, dist: $dist, min at $(argmin(dist)) -> $(ibest+i)")
    mindist = dist[ibest] == 0 ? 0.00001 : dist[ibest]
    jbest = ibest+i
    #disttot = sum((@. haversine_d(mdf.lat[i:jbest-1],mdf.lon[i:jbest-1],mdf.lat[i+1:jbest],mdf.lon[i+1:jbest]) ))#MatoParsing.distance_df(mdf[i:jbest,:]))
    disttot = sum(distance_df(mdf, i, jbest))
    
    tdiff = mdf[jbest,time_col]-mdf[i,time_col]
    ## new speed
    newspeed = _speed(mindist, tdiff) #mindist *3.6 / tdiff
    oldspeed = _speed(disttot, tdiff)

    #println("distlong: $disttot, mindist: $mindist,  ratio: $(disttot/mindist), newspeed $(round(newspeed;digits=1)) oldspeed $(round(oldspeed;digits=1))")
    
    if ( (disttot/mindist) > 15 || (oldspeed > speed_thrs) )&& (newspeed < speed_thrs)  
        ## cut trace
        keep[i+1:jbest-1] .= false
        println(f"Remove {jbest-i-1}-jump of {disttot:4.1f} m -> {mindist:3.1f}, old speed: {oldspeed:3.1f} new speed: {newspeed:3.1f}")
        if tdiff > 60
            println("Warning: removing loop with wrong speed creates time delta of $tdiff s")
        end
    end
end

function remove_large_jumps_df(mdf::AbstractDataFrame, dist::Vector{T}, speed_thrs::Number; time_col::Symbol=:timerec) where T<:AbstractFloat

    #dist = MatoParsing.distance_df(mdf)
    lendf = size(mdf,1)
    if lendf < 3
        return mdf
    end
    times = mdf[2:end,time_col] .- mdf[1:end-1,time_col]
    speed = _speed(dist,times) #@. dist *3.6 / times;
    #dist_thresh = max(500.0,quantile(dist,0.95))
    if quantile(dist,0.95) < 500
        dist_thresh = 500
    else
        dist_thresh = 1000
    end

    mk =@. (speed > speed_thrs) | (dist > dist_thresh) 
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
            disttot = sum(distance_df(mdf, i, jcut))
            tdiff = mdf[jcut,time_col]-mdf[i, time_col]
            avg_speed = _speed(disttot, tdiff)
            if (avg_speed > speed_thrs) && (tdiff < 60)
                keep[i+1:jcut] .= false
                println(f"Remove sec with speed {avg_speed:3.1f} and dist {disttot:3.1f}")
            end
        else
            _check_cut_at_index!(mdf, i, jend, keep, speed_thrs, time_col)
        end
        
    end

    mdf[keep,:]
end

function calc_avg_speed_veh_df(df::AbstractDataFrame)
    idcs_change = findall((df.secIdx[2:end] .-df.secIdx[1:end-1]) .> 0)
    idcs_change = [idcs_change.+1; length(df.secIdx)]
    prev_i =1
    changes = Dict{Symbol,Any}[]
    for ii in idcs_change
        dist = distance_df(df, prev_i, ii)
        totdist = sum(dist)
        tdiff = df.timerec[ii]-df.timerec[prev_i]
        avg_speed = _speed(totdist, tdiff)

        push!(changes,Dict(:section => df.secIdx[prev_i],:distance => totdist, :tottime => tdiff, :speed_avg => avg_speed, :start_t => df.timerec[prev_i],
                :sec_len => ii-prev_i+1))

        prev_i = ii
    end

    DataFrame(changes)
end

function calc_avg_speed_section_df(df::AbstractDataFrame, polyLong::AbstractDataFrame; timecol::Symbol=:timerec)
    idcs_change = findall((df.secIdx[2:end] .-df.secIdx[1:end-1]) .> 0)
    idcs_change = [idcs_change.+1; length(df.secIdx)]
    prev_i =1
    changes = Dict{Symbol,Any}[]
    #pp = df[1,:]
    for ii in idcs_change
        iprev_poly = df[prev_i,:polyPoint]
        inext_poly = df[ii,:polyPoint]
        #println("poly iprev: $iprev_poly inext: $inext_poly")
        dist = distance_df(polyLong, iprev_poly, inext_poly)
        totdist = sum(dist)
        tdiff = df[ii,timecol]-df[prev_i,timecol]
        avg_speed = _speed(totdist, tdiff)

        push!(changes,Dict(:section => df.secIdx[prev_i],:distance => totdist, :tottime => tdiff, :speed_avg => avg_speed, :start_t => df[prev_i,timecol],
                :sec_len => ii-prev_i+1))

        prev_i = ii
    end

    DataFrame(changes)
end

function is_order_df_correct(dfk1::AbstractDataFrame, dfk2::AbstractDataFrame)
    m1=dfk1.timerec[1]
    M1=dfk1.timerec[end]
    m2=dfk2.timerec[1]
    M2=dfk2.timerec[end]
    order = 1
    if (m1 > M2) && (M1 > m2)
        order = 1
    elseif (m2>M1)&&(M2>m1)
        order = -1
    else
        order = 0
    end
    order
end

function add_stats_date(speeddf::AbstractDataFrame; offset_hours::Int=0)
    dates = Dates.unix2datetime.(speeddf.start_t).+Hour(offset_hours)

    speednewdf=hcat(speeddf,
    DataFrame(map(x-> (Dates.monthday(x)..., Dates.hour(x), Dates.minute(x), Dates.dayofweek(x)),dates), [:month,:day,:hour,:minutes,:dayofweek] )
    )
end

"""
Check if the vehicle is going in the opposite direction (from higher section to lower section)
"""
function is_direction_opposite(dftrip)
    diffIdx = dftrip.secIdx[2:end] .- dftrip.secIdx[1:end-1]
    #all(diffIdx.<=0) 
    (sum(diffIdx.<0) > sum(diffIdx .> 0) ) && (dftrip.secIdx[1] >= dftrip.secIdx[end])
end

"""
Divide DataFrame trajectory if new jumps are created
"""
function divide_large_jumps(dftrip::AbstractDataFrame, diff_sec::Integer=2)
    secDiff=dftrip.secIdx[2:end]-dftrip.secIdx[1:end-1]
    icut = findall(secDiff .> diff_sec).+1
    dfs = AbstractDataFrame[]
    if length(icut) > 0
        is = 1
        for i in icut
            push!(dfs,dftrip[is:i-1,:])
            is = i
        end
        push!(dfs,dftrip[is:end,:])
    else
        push!(dfs,dftrip)
    end
    dfs
end