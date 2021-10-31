
curvedLines = function(loc, outDir)
{
    library(data.table)
    library(geosphere)
    
    if(ncol(loc) != 3)
        stop('Input file should have 3 columns named as Lat, Long, name')
    if(!dir.exists(outDir))
        system(sprintf('mkdir %s', outDir))
    
    firstLocations = loc$name
    for(i in 1:length(firstLocations))
    {
        p1_lat = loc$Lat[i];p1_lng = loc$Long[i]
        for(j in 1:length(firstLocations))
        {
            if(j==i)
            {
                p2_lat = p1_lat+0.15;p2_lng = p1_lng+(-0.15)
                a = cbind(p1_lng, p1_lat)
                b = cbind(p2_lng, p2_lat)
                kk = midPoint(a, b)
                r = pointDistance(a, kk, lonlat=TRUE) + pointDistance(b, kk, lonlat=TRUE)
                number_of_ponts = 50
                outFile = paste0(firstLocations[i],'_',firstLocations[j])
                courbature = 4000
                
                #from here run the python code GetEllipsePoints.py with the values for p1_lat, p1_lng, p2_lat, p2_lng,r,number_of_ponts, outDir, outFile parameters
                system(sprintf("python3 GetEllipsePoints.py %s %s %s %s %s %s %s %s %s", p1_lat, p1_lng, p2_lat, p2_lng, r, number_of_ponts, outDir, outFile, courbature))
                
            }
            else
            {
                next
                p2_lat = loc$Lat[j];p2_lng = loc$Long[j]
                a = cbind(p1_lng, p1_lat)
                b = cbind(p2_lng, p2_lat)
                kk = midPoint(a, b)
                r = pointDistance(a, kk, lonlat=TRUE) + pointDistance(b, kk, lonlat=TRUE)
                number_of_ponts = 2000
                outFile = paste0(firstLocations[i],'_',firstLocations[j])
                courbature = 8000

                #from here run the python code GetEllipsePoints.py with the values for p1_lat, p1_lng, p2_lat, p2_lng,r,number_of_ponts, outDir, outFile parameters
                system(sprintf("python3 GetEllipsePoints.py %s %s %s %s %s %s %s %s %s", p1_lat, p1_lng, p2_lat, p2_lng, r, number_of_ponts, outDir, outFile, courbature))
                
            }
        }
    }
}


