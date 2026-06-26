%% ================================================================
%  utils/bootstrap_CI.m
%
%  Compute bootstrap mean and 95% confidence interval for a
%  matrix of observations.
%
%  INPUT
%    X     : [Nevents x Ntime]  matrix
%    Nboot : number of bootstrap resamples
%    alpha : significance level (e.g. 0.05 for 95% CI)
%
%  OUTPUT
%    mu : [1 x Ntime]  column-wise mean
%    lo : [1 x Ntime]  lower CI bound
%    hi : [1 x Ntime]  upper CI bound
%
% ================================================================

function [mu, lo, hi] = bootstrap_CI(X, Nboot, alpha)

    mu = mean(X, 1, 'omitmissing');

    bootMeans = zeros(Nboot, size(X, 2));

    for b = 1:Nboot
        idx = randi(size(X, 1), size(X, 1), 1);
        bootMeans(b, :) = mean(X(idx, :), 1);
    end

    lo = prctile(bootMeans, 100*(alpha/2));
    hi = prctile(bootMeans, 100*(1 - alpha/2));

end
