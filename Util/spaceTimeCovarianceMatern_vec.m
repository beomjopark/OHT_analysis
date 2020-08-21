%% Matern Covariance function 
%% with nu = 3/2 : 1 times differentiable sample path

function cov = spaceTimeCovarianceMatern_vec(lat,long,time,thetas,thetaLat,thetaLong,thetat)
    distSpaceTime = sqrt( ((lat - lat') ./ thetaLat).^2 +...
        ((long - long') ./ thetaLong).^2 +...
        ((time - time') ./ thetat).^2 );

    term1 = sqrt(3).* distSpaceTime;
    cov = thetas .* (1 + term1) .* exp(- term1);
 