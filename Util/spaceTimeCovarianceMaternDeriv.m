function r = spaceTimeCovarianceMaternDeriv(lat1,long1,t1,lat2,long2,t2,thetas,thetaLat,thetaLong,thetat,targetVar)
    %% Derivative of Matern covariance with nu=3/2

    switch targetVar
        case 'lat'
            targetTheta = thetaLat;
            targetDelta = (lat1'-lat2);
        case 'lon'
            targetTheta = thetaLong;
            targetDelta = (long1'-long2);
    end

    distSpaceTime = sqrt(((lat1'-lat2)./thetaLat).^2 + ((long1'-long2)./thetaLong).^2 + ((t1'-t2)./thetat).^2);
    
    r = -3 .* thetas ./ (targetTheta.^2) .* targetDelta .* exp(- sqrt(3).* distSpaceTime);
