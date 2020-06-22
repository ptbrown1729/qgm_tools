function chi = fg_compress_static(beta, mu, nsites)
 % limiting form of fg polarizability/density-density response function
 % I believe this function works much better than the NonIntFG class
 % function for compressibility...
 % In the high temperature limit, expect chi*T = n(1-n) for one component
 % For an m-component gas, this would be chi*T = m*(n/m)(1-n/m) = n*(1-n/m)
 
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
    nsites = 100;
end

if numel(beta) ~= 1
% if we have multiple (beta, mu) points to evaluate, loop through each and
% call the function on those individual.
% TODO: probably faster to do this vectorized, but need to recognize when
% we are creting very large arrays.
    chi = zeros(size(beta));
    for ii = 1 : numel(beta)
        chi(ii) = fg_density( beta(ii), mu(ii), nsites);
    end
else
    % allowed momenta
%     ks = (2*pi) / nsites * (0 : nsites -1 );
%     [kys, kxs] = ndgrid(ks, ks);
    [kxs, kys] = get_allowed_kvects(nsites, nsites, 'balanced');
    % dispersion
    Es = -2 * (cos(kxs) + cos(kys));

    % fermi function
    fermi = @(b, mu, e) 1 ./ (exp(b * (e - mu)) + 1);

    chi = ...
        beta * 1/nsites^2 * squeeze(sum(sum(...
        fermi(beta, mu, Es).^2 .* exp(beta * (Es - mu)), 2), 1));
end

end