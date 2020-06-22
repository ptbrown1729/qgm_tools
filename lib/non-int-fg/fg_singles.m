function [n_singles] = fg_singles(beta, mu_up, mu_dn, nsites)
% compute the singles density of a non-interacting fermi gas with a tightbinding
% dispersion at given value of beta and mu. All energies are 
% measured in fractions of the hopping energy t.
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

if ~exist('nsites', 'var')
    nsites = 50;
end

% computation
n_ups = fg_density(beta, mu_up, nsites);
n_dns = fg_density(beta, mu_dn, nsites);
n_singles = n_ups + n_dns - 2 * n_ups .* n_dns;

end

