function crps = crps_normal(mean, sd, pred)
    z = (pred - mean) ./ sd;
    crps = sd .* ((z .* (2 * normcdf(z) - 1)) + 2 * normpdf(z) - 1/ sqrt(pi));
end