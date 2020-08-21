function cov = spaceTimeCovarianceMaternDeriv_vec(lat,long,time,thetas,thetaLat,thetaLong,thetat,targetVar)
%% Derivative of Matern Covariance function 
%% with nu = 3/2 : 1 times differentiable sample path

    switch targetVar
    case 'lat'
        targetTheta = thetaLat;     
        targetTerm = ((lat - lat') ./ thetaLat).^2;
    case 'lon'
        targetTheta = thetaLong;
        targetTerm = ((long - long') ./ thetaLong).^2;
    end

    distSpaceTime = sqrt( ((lat - lat') ./ thetaLat).^2 +...
        ((long - long') ./ thetaLong).^2 +...
        ((time - time') ./ thetat).^2 );
    cov = 3 .* thetas ./ (targetTheta.^2) .* (1  - sqrt(3) .* targetTerm ./ distSpaceTime) .* exp(- sqrt(3) .* distSpaceTime);
    cov(1:(numel(lat)+1):end) = 3 .* thetas ./ (targetTheta.^2);
