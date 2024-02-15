struct Stop{T<:Real}
    lat::T
    lon::T
    name::String
end

function Stop(dd::Dict{String,Any},namekey::String="gtfsId")
    Stop(dd["lat"],dd["lon"],dd[namekey])
end