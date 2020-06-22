function [dens] = fg_density(beta, mu, nsites)
% compute the density of a non-interacting fermi gas with a tightbinding
% dispersion at given value of beta and mu. All energies are 
% measured in fractions of the hopping energy t.
%
% beta and mu may be arrays of any size or number of dimensions, but they
% should either have the same size, or one of them may be an arbitrary array
% and the other a single number.
%
% nsites is the number of sites in each direction to include in the
% problem.

if numel(beta) == 1 && numel(mu) > 1
    beta = beta * ones(size(mu));
end

if numel(beta) > 1 && numel(mu) == 1
    mu = mu * ones(size(beta));
end

if ~isequal( size(beta), size(mu))
    error('beta and mu must be the same size');
end

if ~exist('nsites', 'var')
    nsites = 50;
end

%

if numel(beta) ~= 1
% if we have multiple (beta, mu) points to evaluate, loop through each and
% call the function on those individual.
% TODO: probably faster to do this vectorized, but need to recognize when
% we are creting very large arrays.
    dens = zeros(size(beta));
    for ii = 1 : numel(beta)
        dens(ii) = fg_density( beta(ii), mu(ii), nsites);
    end
else
    % allowed k-vectors
%     ks = (2*pi) / nsites * (0 : nsites - 1 );
%     [kys, kxs] = ndgrid(ks, ks);
    [kxs, kys] = get_allowed_kvects(nsites, nsites, 'balanced');
    % lattice dispersion 
    Es = -2 * (cos(kxs) + cos(kys));
    Es = Es(:);

    % fermi function
    fermi = @(b, mu, e) 1 ./ (exp(b * (e - mu)) + 1);

    % total size of region is nsites^2, so must divide by this to get
    % density
    dens = 1 / nsites^2 * squeeze( sum( fermi(beta, mu, Es) ) );
end

end

