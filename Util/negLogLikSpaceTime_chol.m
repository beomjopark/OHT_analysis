function val = negLogLikSpaceTime_chol(params,profLatAggr,profLongAggr,profJulDayAggr,intResidAggr,kernelType)

thetas = exp(params(1));
thetaLat = exp(params(2));
thetaLong = exp(params(3));
thetat = exp(params(4));
sigma = exp(params(5));

%disp([thetas,thetaLat,thetaLong,thetat,sigma]);

largeVal = 1e10;

if ~(all([thetas,thetaLat,thetaLong,thetat,sigma] < largeVal))
    val = NaN;
    disp(val);
    return;
end

if ~(all([thetas,thetaLat,thetaLong,thetat,sigma] > eps))
    val = NaN;
    disp(val);
    return;
end

nYear = size(intResidAggr,2);

val = 0;

for iYear = 1:nYear
    %disp(iYear);
    
    profLatYear = profLatAggr{iYear};
    profLongYear = profLongAggr{iYear};
    profJulDayYear = profJulDayAggr{iYear};
    intTempResYear = intResidAggr{iYear};
    
    nRes = length(intTempResYear);
    if ~nRes
        continue;
    end
    %disp(nRes);
    
    covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                profLatYear, profLongYear, profJulDayYear, thetas,thetaLat,thetaLong,thetat);
%{
    covObs = zeros(nRes,nRes);
    
    %tic;
    for i = 1:nRes
        %disp(i);
        for j = 1:nRes
            covObs(i,j) = spaceTimeCovarianceExpGeom(profLatYear(i),profLongYear(i),profJulDayYear(i),profLatYear(j),profLongYear(j),profJulDayYear(j),thetas,thetaLat,thetaLong,thetat);
        end
    end
    %toc;
%}   
    %tic;
    K = covObs + sigma.^2*eye(nRes);
    K = (K + K') ./ 2;
    val = val + 2*sum(log(diag(chol(K)))) + (intTempResYear)'*(K\intTempResYear) + log(2*pi)*nRes;
    %toc;
    
end

val = 0.5*val;

%disp(val);