module MatoParsing

using DataFrames
using PyFormattedStrings
import Base: length

export haversine_d, length

length(x::AbstractDataFrame) = size(x,1)

include("haversine.jl")
include("trace_parsing.jl")

end # module MatoParsing
