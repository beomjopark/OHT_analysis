function r = spaceTimeCovarianceMatern(lat1,long1,t1,lat2,long2,t2,thetas,thetaLat,thetaLong,thetat)

distSpaceTime = sqrt(((lat1'-lat2)./thetaLat).^2 + ((long1'-long2)./thetaLong).^2 + ((t1'-t2)./thetat).^2);

term1 = sqrt(3).* distSpaceTime;
r = thetas .* (1 + term1) .* exp(- term1);