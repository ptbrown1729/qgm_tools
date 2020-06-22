function [density, phase_space_density] = getPhaseSpaceDensity(N, T, sr, sz, Constants)
% compute phase space density assuming a gaussian density profile

%Get constants
mLi = Constants.mass;
K = Constants.K;
h = Constants.h;


% N=800e6;
% T=40e-6;
magnification = 30;
pix_size = 6.5e-6;
nbin = 8;

% sigmas
sigmar = sr * nbin * pix_size / magnification;
sigmaz = sz * nbin * pix_size / magnification;
% TODO: why is this here? Seems like this function is broken?
sx = 5e-6;
% atom density
volume_gaussian = sigmaz * sigmar * sx * (2*pi)^(3/2);
density = N / volume_gaussian;
% phase space density
lambda_db = h / sqrt(2 * pi * mLi * K * T);
phase_space_density = density * lambda_db ^ 3;
