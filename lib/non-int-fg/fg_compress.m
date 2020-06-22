function [kx, ky, chi] = fg_compress(omega, beta, mu, nsites)
% compute charge susceptibility (polarizability) in a free Fermi gas.
% This is the retarded density-density correlation function in momentum
% and frequency space.
% G^R(\rho(q), \rho(-q)) = -i/hbar theta(t) <[rho(q, t), rho(-q, 0)]>
% where \rho(q) = \sum_k c^\dag_k c_{k+q}, which is the fourier
% transform of n_i = c^\dag_i c_i
% G(q, \omega) is the Lindehard function, as can be shown by computing
% the Matsubara green's function, applying Wick's theorem, and reducing.
% See for example "Many-body quantum theory in condensed matter physics"
% by Henrik Bruus and Karsten Flensberg, section 10.7
   
if ~isequal( size(beta), size(mu))
    error('beta and mu must be the same size');
end
   
if ~exist('nsites', 'var')
    nsites = 100;
end

% allowed momenta 
ks = (2*pi) / nsites * (0 : nsites -1 );
[kys, kxs, k_shift_y, k_shift_x] = ndgrid(ks, ks, ks, ks);

% dispersions
Es = -2 * (cos(kxs) + cos(kys));
Es_shifted = -2 * (cos(kxs - k_shift_x) + cos(kys - k_shift_y));

% fermi function
fermi = @(b, mu, e) 1 ./ (exp(b * (e - mu)) + 1);

denominator = (omega + Es - Es_shifted);
% if this is negative, then the excitation is not allowed
denominator(denominator < 0) = NaN;

chi = ...
1 / beta * 1/nsites^2 * squeeze(sum(sum(...
   (fermi(beta, mu, Es) - fermi(beta, mu, Es_shifted)) ./ ...
   denominator, 2), 1));

kx = squeeze(k_shift_x(1, 1, :, :));
ky = squeeze(k_shift_y(1, 1, :, :));
    
end
 
