module MatoParsing

using DataFrames
using PyFormattedStrings
using StatsBase
import Base: length

export haversine_d, length, distance_df, linspace

length(x::AbstractDataFrame) = size(x,1)
function linspace(s::Number,e::Number,np::Integer)
    spacing = (e-s)/(np-1)
    collect(s:spacing:e)
end
include("io.jl")

include("haversine.jl")
include("trace_parsing.jl")
include("polyline.jl")

end # module MatoParsing
