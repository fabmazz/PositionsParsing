module MatoParsing

using DataFrames
using PyFormattedStrings
using StatsBase
import Base: length

export haversine_d, length, distance_df

length(x::AbstractDataFrame) = size(x,1)
include("io.jl")

include("haversine.jl")
include("trace_parsing.jl")
include("polyline.jl")

end # module MatoParsing
