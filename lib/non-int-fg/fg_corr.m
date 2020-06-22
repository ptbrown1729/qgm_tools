function [corr] = fg_corr(beta, mu, index, nsites)
% compute the density correlations = <n_i n_j> - <n_i><n_j> of a 
% non-interacting fermi gas with a tightbinding
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

if ~exist('index', 'var')
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
    corr = zeros(size(beta));
    for ii = 1 : numel(beta)
        corr(ii) = fg_corr( beta(ii), mu(ii), index, nsites);
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

    ft = exp(1i * kxs * index(1) + 1i * kys * index(2));
    %now compute density/site by summing over allowed states
    corr = -abs(1/nsites^2 * sum(sum( ft .* fermi(beta, mu, Es) ))).^2;
end

end

