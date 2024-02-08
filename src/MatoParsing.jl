module MatoParsing

using DataFrames
using PyFormattedStrings
using StatsBase
import Base: length
using Dates

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
include("sections_stops.jl")

include("processing.jl")

end # module MatoParsing
