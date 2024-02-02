
# Define the Haversine distance function
function haversine_d(lat1::T, lon1::T, lat2::T, lon2::T, radius::AbstractFloat=6_371_000.0) where T<: AbstractFloat
    # Earth radius in kilometers
    #const R::Float64 = 6_371_000
    
    # Convert latitude and longitude from degrees to radians
    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)
    
    # Haversine formula
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    
    # Distance in kilometers
    distance::T = radius * c
    return distance
end

