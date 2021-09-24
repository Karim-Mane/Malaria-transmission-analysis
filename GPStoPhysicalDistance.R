
#----- This function uses the GPS coordinates of two location and returns the physical distance between them.

convertGPSintoPhysicalDistance = function(long1, lat1, long2, lat2)  
{
    deltaLong=deltaLat=a=b=b2=c=d=var1=var2=kDistance=0
    p = pi/180
    deltaLong = (long2-long1)*p
    deltaLat = (lat2-lat1)*p
    a = sin(deltaLat/2.0)
    b = sin(deltaLong/2.0)
    b2 = b^2
    c = cos(lat1*p)
    d = cos(lat2*p)
    var1 = a^2 + c*d*b2
    var2 = 2*atan2(sqrt(var1), sqrt(1-var1))
    kDistance = round(6370*var2, digits = 0)
    return(kDistance)
}