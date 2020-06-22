function [doubles_corr] = fg_doubles_corr(beta, mu_up, mu_dn, index, nsites)
% compute the singles correlations of a non-interacting fermi gas with a tightbinding
% dispersion at given value of beta and mu. All energies are 
% measured in fractions of the hopping energy t.
%
% The singles correlator can be calculated from density correlations using
% Wick's theorem according to
% <d(i) d(j)>_c = n_up^2 * <n_up(i) n_up(j)>_c 
%                 n_dn^2 * <n_dn(i) n_dn(j)>_c 
%                 + <n_up(i) n_up(j)>_c * <n_dn(i) n_dn(j)>_c 
%               
% beta and mu may be arrays of any size or number of dimensions, but they
% should either have the same size, or one of them may be an arbitrary array
% and the other a single number.
%
% nsites is the number of sites in each direction to include in the
% problem.

if numel(beta) == 1
    if numel(mu_up) > 1
        beta = beta * ones(size(mu_up));
    elseif numel(mu_dn) > 1
        beta = beta * ones(size(mu_up));
    end
end

if numel(mu_up) == 1 
    if numel(beta) > 1
        mu_up = mu_up * ones(size(beta));
    elseif numel(mu_dn) > 1
        mu_up = mu_up * ones(size(mu_dn));
    end  
end

if numel(mu_dn) == 1
    if numel(beta) > 1
        mu_dn = mu_dn * ones(size(beta));
    elseif numel(mu_up) > 1
        mu_dn = mu_dn * ones(size(mu_up));
    end
end

if ~isequal( size(beta), size(mu_up)) || ~isequal( size(beta), size(mu_dn))
    % this also force size(mu_up) == size(mu_dn)
    error('beta and mu must be the same size');
end

if ~exist('index', 'var') || isempty(index)
    index = [0, 1];
end

if ~exist('nsites', 'var')
    nsites = 50;
end

if numel(beta) ~= 1
% if we have multiple (beta, mu) points to evaluate, loop through each and
% call the function on those individual.
% TODO: probably faster to do this vectorized, but need to recognize when
% we are creting very large arrays.
    doubles_corr = zeros(size(beta));
    for ii = 1 : numel(beta)
        doubles_corr(ii) = fg_doubles_corr( beta(ii), mu_up(ii), mu_dn(ii), index, nsites);
    end
else
    % allowed k-vectors
%     ks = (2*pi) / nsites * (0 : nsites - 1 );
%     [kys, kxs] = ndgrid(ks, ks);
    [kxs, kys] = get_allowed_kvects(nsites, nsites, 'balanced');
    % lattice dispersion 
    Es = -2 * (cos(kxs) + cos(kys));


    % fermi function
    fermi = @(b, mu, e) 1 ./ (exp(b * (e - mu)) + 1);

    % now compute density/site by summing over allowed states
    n_up = 1 / nsites^2 * squeeze( sum( sum(fermi(beta, mu_up, Es), 2), 1));
    n_dn = 1 / nsites^2 * squeeze( sum( sum(fermi(beta, mu_dn, Es), 2), 1));
    
    % compute correlators
    ft = exp(1i * kxs * index(1) + 1i * kys * index(2));
    corr_up = -abs(1/nsites^2 * sum(sum( ft .* fermi(beta, mu_up, Es) ))).^2;
    corr_dn = -abs(1/nsites^2 * sum(sum( ft .* fermi(beta, mu_dn, Es) ))).^2;
    
    % combine using Wick's theorem
%     doubles_corr = n_up^2 * corr_up + n_dn^2 * corr_dn + corr_up * corr_dn;
    doubles_corr = n_up^2 * corr_dn + n_dn^2 * corr_up + corr_up * corr_dn;
end

end

