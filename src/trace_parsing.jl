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


function row_avg(df, j, i, count)
    for key in [:lat, :lon, :heading]
        if (k=="heading") && ismissing(df[j,key])
            df[j,k] = df[i,k] ## ignore None heading
        else
            df[j,k] = (df[j,k]*count[j] + df[i,k])/(count[j]+1)
        end
    end
    count[j]+=1
    count[i]=0
end

function average_pos_equaltime(spec_df)
    lendf = size(spec_df,1)
    count = ones(Int,lendf) #np.ones(len(spec_df),dtype=np.int_)
    secs_diff = spec_df.timerec[2:end].-spec_df.timerec[1:end-1];

    idcs_eq = findall(secs_diff==0) #np.where(secs_diff == 0)[0]
    if length(idcs_eq) == 0
        return spec_df
    else
        println("Have $(length(idcs_eq)) value with same time")
    end
    for i in idcs_eq
        for j=(i-1):-1:1
            if count[j] > 0
                row_avg(spec_df, j, int(i), count)
                break
            end
        end
    end
    
    return spec_df[count.>0,:]
end