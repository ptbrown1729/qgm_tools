function mu = fg_mu(beta, n, nsites)
% compute the chemical potential of a non-interacting fermi gas with a tightbinding
% dispersion at given value of beta and n. All energies are 
% measured in fractions of the hopping energy t.
%
% beta and n may be arrays of any size or number of dimensions, but they
% should either have the same size, or one of them may be an arbitrary array
% and the other a single number.
%
% nsites is the number of sites in each direction to include in the
% problem.
if numel(beta) == 1 && numel(n) > 1
    beta = beta * ones(size(n));
end

if numel(beta) > 1 && numel(n) == 1
    n = n * ones(size(beta));
end

if ~isequal( size(beta), size(n))
    error('beta and n must be the same size');
end

if any(n > 1) || any(n < 0)
    error('n must be between 0 and 1');
end

if ~exist('nsites', 'var')
    nsites = 100;
end

if numel(beta) ~= 1
% if we have multiple (beta, mu) points to evaluate, loop through each and
% call the function on those individual.
% TODO: probably faster to do this vectorized, but need to recognize when
% we are creting very large arrays.
    mu = zeros(size(beta));
    for ii = 1 : numel(beta)
        mu(ii) = fg_mu( beta(ii), n(ii), nsites);
    end
    
else
    fn = @(mu) fg_density(beta, mu, nsites) - n;

    init_guess = 0;

    try
        if n > 0.01
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective', 'display', 'off');
        else
            % at very low densities it seems like the levenberg-marquardt
            % algorithm converges better, but it is also slower.
            options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt', 'display', 'off');
        end
        [mu, resnorm, residual, exitflat, output, lambda, jacobian] = lsqnonlin(fn, init_guess, [], [], options);
    catch
        mu = NaN;
    end
        
    % ninety_five_per_confint = nlparci(mu, residual, 'jacobian', jacobian);  
    % mu_unc = (ninety_five_per_confint(:,2)-ninety_five_per_confint(:,1))/3.92;
end

end

