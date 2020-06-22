function [tTB, tTBx, tTBy, tTBdiag, u,...
    band_widths, GapFromGroundBand, BandTopToGroundBandBottom, ...
    MeanSplittingToGroundBand, EBand, ...
    WannierFn, latticeAngle] = ...
    Lattice_2D_Asymmetric(s, as, zConf, Atten, latticeAngle, WannierCalcOn, show_plots)
% [tTB, tTBx, tTBy, tTBdiag, u,...
%     BandWidths, GapFromGroundBand, BandTopToGroundBandBottom, ...
%     MeanSplittingToGroundBand, EBand, ...
%     WannierFn, latticeAngle] = ...
%     Lattice_2D_Asymmetric(s,as,zConf,Atten,WannierCalcOn,ShowPlots)
%
% Solve lattice band structure matrix problem, which is
% hbar^2/(2m) * (k - \Kappa)^2 c_{k-\Kappa} + \sum_K' U(K'-\Kappa) c_{k-K'} = \epsilon c_{k-\Kappa}
% where \Kappa are the reciprocal lattice vectors, and k is a vector in the
% Brillouin zone.
%
% We take this as an eigenvalue problem, regarding the first term as a
% diagonal matrix indexed by \Kappa (for fixed k in the Brillouin zone).
% We define a matrix U_{K', \Kappa} = U(K' - \Kappa), where U is the
% spatial fourier transform of the periodic potential.
%
% We compute the tunneling element from the band structure assuming the
% tight-binding limit.
% 
% We compute the interaction energy U from the wannier function, the
% vertical confinement, and the scattering length
%
%%% Arguments
% s: lattice depth in Ers
%
% as: scattering length in Bohr radii
%
% zConf: is axial trapping frequency in Hz
%
% Atten: is factor of attenuation for lattice retro = 0.464 or so. 
%
% WannierCalcOn: is a boolean which specifies whether or not to calculate
% Wannier functions.
%
% ShowPlots: is a boolean which specified whether or not to display summary
% plots.
%
%%% Return values
% return values are given in units of Er
%
% tTB: is 1/8 the bandwidth of the ground band.
%
% tTBx: x direction tunneling energy in Ers
%
% tTBy:
%
% tTBdiag:
%
% u: is the interaction energy
%
% Bandwidths: the widths of each band.
%
% GapFromGroundBand: is the distance between the top of the ground band and
% the bottom of the given band
%
% BandTopToGroundBandBottom: is the distance between the bottom of the ground
% band and the top of the given band
%
% MeanSplittingToGroundBand: is the mean distance between the ground band and
% the given band (also the distance between the average energies of the
% band)
%
% Eband: is an Nbands x Nkpts x Nkpts array giving the energy of the bands in
% the Brillouin zone.
%
% latticeAngle: is the angle between the two lattice beams. This can be
% determined from the reconstruction files. It effects the spacing of the
% two axes of the lattice, but generally does not change the angle between the lattice axes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set number of reciprocal lattice vectors, Brillouin zone vectors, etc. to
% include in calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of reciprocal lattice vectors.
% i.e. Fourier components of U. Increase for more accurate lattice.
n = 10;
nn_recp_vects = (2 * n + 1)^2;

% number of Brillouin zone vectors
qsteps = 12;
nn_bz_vects = 2 * qsteps + 1;
i_center = qsteps + 1; % center index

% number of bands to keep
bands_to_store = 8;

% wannier function real space grid settings
real_space_steps = 10; % number of points per unit cell
nWZ = 3;  % number of WZ cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('WannierCalcOn', 'var') || isempty(WannierCalcOn)
    WannierCalcOn = 1;
end

if ~exist('Atten', 'var') || isempty(Atten)
    Atten = 1;
end

if ~exist('latticeAngle', 'var') || isempty(latticeAngle)
    % theta = 2 * atan(ax /ay) where ax and ay are the lattice spacings in the two directions
    latticeAngle = 91.6267; % angle btw lattice beams in deg.  
end
lattice_spacing_ratio = tan(latticeAngle / 2 * pi/180);

abohr = 5.29e-11;
if ~exist('as', 'var') || isempty(as)
    asc = 1000 * abohr;
else
    asc = as * abohr;
end

if ~exist('zConf', 'var') || isempty(zConf)
    zConf = 20e3;
end

if ~exist('show_plots', 'var') || isempty(show_plots)
    show_plots = 1;
end

%define useful constants
hbar = 1.055e-34; 
m = 6.015*1.66e-27;
a = 1064e-9 / (2^0.5); 
Er = (pi * hbar)^2 / (2 * m * a^2);

% alpha = cos(latticeAngle / 90 * (pi/2) / 2);
% beta = sin(latticeAngle / 90 * (pi/2) / 2);
alpha = cos(0.5 * latticeAngle * pi/180);
beta = sin(0.5 * latticeAngle * pi/ 180);

% Expression for real space and arbitrary polarization
% Here \alpha is the angle of the polarization vector wrt to vertical direction
% The beams come at directions (x+y) and (x-y), and the lattice principle
% axes are long x and y.
% V1 = s/4 * cos(\alpha)^2; 
% V2 = -0.5 * abs(V1)
% V(x,y) = -V1[cos(kx)+cos(ky)] + V2[cox(kx+ky)+cos(kx-ky)]
%
% Expression along principle axes of lattice, supposing vertical
% polarization
% V(x,y) = V0 * [ 1 - (1 + r^2 + 2*r*cos( 2*k*x*cos(theta/2) ) ) / (1 + 2r
% + r^2) * 0.5 * ( 1 + cos( 2*k*y*sin( theta/2 ) ) )] 
% k = 2pi/a
% at least before accounting for beam angles...

% reciprocal lattice vectors
% build ixy to index c's by G value. G's run from -n to n
recp_latt_vects = transpose(-n : n);

% get all reciprocal lattice vectors as a 1D list
% TODO: mixed up x and y's intentionally to match what was previously
% written. Should correct this in the end though. Maybe entire script is
% consistent this it right now?
[iyiy, ixix] = meshgrid(recp_latt_vects, recp_latt_vects);
ix = ixix(:);
iy = iyiy(:);

% Now get differences between reciprocal lattice vectors
% dkx_ij = ix(j) - ix(i)
[aa, bb] = meshgrid(ix, ix);
dQx = bb - aa;

% dky_ij = iy(j) - ix(i)
[aa, bb] = meshgrid(iy, iy);
dQy = bb - aa;

% Fourier components of the potential

% Construct the potential energy part of the Hamiltonian, U
% i.e. the part of the Hamiltonian that is momentum independent
% Lattice fourier components.
% depthscale = (Atten+(Atten+1)^2/4)/2; % adjusted lattice depth (due to attenuation)
V1 = (3/2) * s / (1 + Atten^2 + 4 * Atten); 
V2 = -1/2 * abs(V1);
V10 = -V1/2 * Atten;
V01 = -V1/2 * (1 + Atten^2) / 2;
V11 = V2/2 * Atten;

% I think these are the right values though...
% V10 = s * 2 * Atten / (4 * (1 + Atten^2 + 2 * Atten) );
% V01 = s * (1 + Atten^2) / (4 * (1 + Atten^2 + 2 * Atten) );
% V11 = - s * Atten / (4 * (1 + Atten^2 + 2 * Atten) );

U = V10 * (abs(dQx) == 1 & dQy == 0) ...
  + V01 * (dQx == 0 & abs(dQy) == 1)...
  + V11 * (abs(dQx) == 1 & abs(dQy) == 1);

% construct q-dependent part of Hamiltonian and diagonalize eigenvalue 
% problem for each q-vector in the BZ
kxs = (-qsteps : qsteps) / qsteps;
kys = (-qsteps : qsteps) / qsteps;
[kxkx, kyky] = meshgrid(kxs, kys);

EBand = zeros(bands_to_store, nn_bz_vects, nn_bz_vects); 
FTpsi = zeros(nn_recp_vects, nn_bz_vects, nn_bz_vects); 
for ii = 1 : length(kxs)
    for jj = 1 : length(kys)
        % construct full hamiltonian from the potential part we calculated above
        % and the kinetic term along diagonal
        H = U + diag( 2 * alpha^2 * (kxs(ii) - 2 * ix).^2 + ...
                      2 * beta^2  * (kys(jj) - 2 * iy).^2);
                  
        % diagonalize problem at this q-vector in the BZ
        [eig_vects, eig_es] = eig(H);
        [EnD, I] = sort(diag(eig_es));
        
        % Record FT of ground band Bloch wavefunction
        % pick out eigenvector w/ lowest E eigenvalue    
        EBand(:, jj, ii) = EnD(1 : bands_to_store);
        FTpsi(:, jj, ii) = eig_vects(:, I(1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get hopping from band structure assuming tightbinding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we suppose band structure is
% E(kx,ky) = -2 [ tx* cos(kx) + ty * cos(ky) + td * cos(kx+ky) + td * cos(kx-ky) ]
% t = 1/8 * bandwidth
%   = 1/8 * [E(pi,pi) - E(0,0)]
% tx = 1/4 * [E(pi,0) - E(0,0)] (no diagonal tunneling)
% tx = 1/8 * ( [E(pi,0) - E(0,0)] + [E(pi,pi) - E(0,pi)]) (finite diagonal tunneling)
% ty = 1/4 * [E(pi,0) - E(0,0)]
% ty = 1/4 * ([E(0, pi) - E(0,0)] + [E(pi,pi) - E(pi,0)])
% td = 1/16 * [E(0, pi) + E(pi, 0) - E(0, 0) - E(pi, pi)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gband = squeeze(EBand(1, :, :));
tTB = 1/8 * ( max(gband(:)) - min(gband(:)) );
tTBx = 1/8 * ( gband(end, i_center) - gband(i_center, i_center) + ...
               gband(end, end) - gband(i_center, end) );
tTBy = 1/8 * ( gband(i_center, end) - gband(i_center, i_center) + ...
               gband(end, end) - gband(end, i_center) );
tTBdiag = 1/16 * ( gband(i_center, end) + gband(end, i_center) - ...
                   gband(i_center, i_center) - gband(end, end) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandwidths and splittings for all bands.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band_widths = max(max(EBand, [], 3), [], 2) - min(min(EBand, [], 3), [], 2); 
GapFromGroundBand = min(min(EBand(2:end, :, :), [], 3), [], 2) - max(max(EBand(1, :, :), [], 3), [], 2);
BandTopToGroundBandBottom = max(max(EBand(2:end, :, :), [], 3), [], 2) - min(min(EBand(1, :, :), [], 3), [], 2);
MeanSplittingToGroundBand = mean(mean(EBand(2:end, :, :), 3), 2) - mean(mean(EBand(1, :, :), 3), 2);

d_band_gaps = [mean(mean(EBand(4, :, :) - EBand(1, :, :))),...
               mean(mean(EBand(5, :, :) - EBand(1, :, :))),...
               mean(mean(EBand(6, :, :) - EBand(1, :, :)))];

Gap0toPMin = GapFromGroundBand(1); 
Gap0toPMax = BandTopToGroundBandBottom(1); 
NextHigherBandsMin = min(min(EBand(7, :, :))) - max(max(EBand(1, :, :)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate Wannier function to get interaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = 0;
WannierFn = 0;
if WannierCalcOn

    % Fix gauge by measuring the euclidean distance between neigboring vectors
    % and choosing a vector sign to minimize it.
    % first fix gauge going along edge of Brillouin zone in the x-direction
    for jj = 1 : length(kxs) - 1
        if sum((FTpsi(:, 1, jj) .* FTpsi(:, 1, jj + 1)), 1) < 0
            FTpsi(:, 1, jj + 1) = -FTpsi(:, 1,jj + 1);
        end
    end
    
    % now fix all the other gauges
    for ii = 1 : length(kys) - 1
        for jj = 1 : length(kxs)
            if sum((FTpsi(:, ii, jj) .*  FTpsi(:, ii + 1, jj)), 1) < 0
                FTpsi(:, ii + 1, jj) = -FTpsi(:, ii + 1, jj);
            end
        end
    end
    
    % Evaluate and sum Bloch wavefunctions over a quarter WZ cell   
    w = zeros(nWZ * real_space_steps, nWZ * real_space_steps); 
    [rr, ss] = meshgrid(kxkx(:), ix);
    [tt, uu] = meshgrid(kyky(:), iy);
    FTpsi_reshaped = reshape(FTpsi, [nn_recp_vects, numel(kxkx)]);
    for jj = 1 : real_space_steps * nWZ
        for ii = 1 : real_space_steps * nWZ
            exp_factor = exp(1i * (jj - 0.5) * pi / real_space_steps * (rr - 2 * ss) + ...
                             1i * (ii - 0.5) * pi / real_space_steps * (tt - 2 * uu) );
            w(ii, jj) = sum(exp_factor(:) .* FTpsi_reshaped(:));
        end
    end

    %calculated for one quadrant. From symmetry, generate full wannier fn.
    wfull= cat(2, flip(cat(1, flip(w, 1), w), 2), cat(1, flip(w, 1), w));
    % normalize
    norm = (sum(sum(conj(wfull) .* wfull, 1), 2));   
    wnorm = wfull / (norm^0.5);
    WannierFn = wnorm;
    
    % Calculate interaction energy in Er
    u0 = 8 * a^2 * asc / pi; %(4*pi*as*hbar^2/m)/Erec
    %First term is harmonic oscillator function in z direction
    u = (m * zConf * 2 * pi / (2 * pi * hbar))^0.5 * u0 *...
        sum(sum((wnorm .* conj(wnorm)).^2, 2), 1) / ((a / (real_space_steps))^2);
    
    %can also check hopping?
%     Vfn = @(x,y) V1*[cos(2*pi*x) + cos(2*pi*y)] + V2*[cos(2*pi*(x+y)) + cos(2*pi*(x-y))];
%     [xx,yy] = meshgrid([-steps*nWZ:-1,1:steps*nWZ],[-steps*nWZ:-1,1:steps*nWZ]);
%     xx = xx/steps; yy = yy/steps;
%     V = Vfn(xx,yy);
%     HoppingVInt = (m*zConf*2*pi/(2*pi*hbar))^0.5*sum(sum(wnorm(:,steps+1:end).*V(:,steps+1:end).*conj(wnorm(:,1:end-steps))))/((a/(steps))^2);
    
%Density dependent hopping corrections
    txDens = (m * zConf * 2 * pi / (2 * pi * hbar))^0.5 * ...
        u0*sum(sum(wnorm(:, real_space_steps+1:end).^2 .* ...
        conj(wnorm(:, real_space_steps+1:end)) .* ...
        conj(wnorm(:, 1:end-real_space_steps)), 2), 1) / ((a / (real_space_steps))^2);
    tyDens = (m * zConf * 2 * pi / (2 * pi * hbar))^0.5 * ...
        u0 * sum(sum(wnorm(real_space_steps+1:end, :).^2 .* ...
        conj(wnorm(real_space_steps+1:end, :)) .* ...
        conj(wnorm(1:end-real_space_steps, :)), 2), 1) / ((a / (real_space_steps))^2);
    tdiagDens = (m * zConf * 2 * pi / (2 * pi * hbar))^0.5 * ...
        u0 * sum(sum(wnorm(real_space_steps+1:end, real_space_steps+1:end).^2 .* ...
        conj(wnorm(real_space_steps+1:end, real_space_steps+1:end)) .* ...
        conj(wnorm(1:end-real_space_steps,1:end-real_space_steps)), 2), 1) / ((a / (real_space_steps))^2);
    
else
    txDens = 0;
    tyDens = 0;
    tdiagDens = 0;
    u = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Print and plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_results = 1;
if print_results
    %print information about lowest bands
    fprintf('#################################\n');
    fprintf('Tight-binding results\n');
    fprintf('#################################\n');
    fprintf('Er = %0.3f KHz \n', Er/(hbar*2*pi)/1e3);
    fprintf('tTB = %0.3f Er = %0.1f Hz\n', tTB, tTB * Er/(hbar*2*pi));
    fprintf('tTBx = %0.3f Er = %0.1f Hz\n', tTBx, tTBx*Er/(hbar*2*pi));
    fprintf('tTBy = %0.3f Er = %0.1f Hz\n', tTBy, tTBy*Er/(hbar*2*pi));
    fprintf('tTBdiag = %0.3f Er = %0.1f Hz\n', tTBdiag, tTBdiag*Er/(hbar*2*pi));

    fprintf('#################################\n');
    fprintf('Wannier function results\n');
    fprintf('#################################\n');
    fprintf('U = %0.3f Er = %0.2f KHz \n', u, u*Er/(hbar*2*pi)/1e3);
    fprintf('U/t = %0.2f \n', u/tTB);
    fprintf('txDens = %0.3f Er = %0.1f Hz\n', txDens, txDens * Er / (hbar*2*pi))
    fprintf('tyDens = %0.3f Er = %0.1f Hz\n', tyDens, tyDens * Er / (hbar*2*pi))
    fprintf('tDens mean/ t mean = %0.3f\n', (txDens + tyDens)/(tTBx+tTBy));

    fprintf('#################################\n');
    fprintf('Lattice geometry information\n');
    fprintf('#################################\n');
    fprintf('Angle between incoming lattice beams = %0.3f deg\n', latticeAngle);
    fprintf('ax = %0.2f lambda \n', alpha);
    fprintf('ay = %0.2f lambda \n', beta);
    fprintf('lattice spacing ratio = %0.2f\n', lattice_spacing_ratio);

    fprintf('#################################\n');
    fprintf('Band splitting information\n');
    fprintf('#################################\n');
    fprintf('Ground to p bands min gap = %0.3f Er = %0.1f KHz \n', Gap0toPMin, Gap0toPMin * Er / (hbar*2*pi) / 1e3);
    fprintf('Ground to p bands max gap = %0.3f Er = %0.1f KHz \n', Gap0toPMax, Gap0toPMax * Er / (hbar*2*pi) / 1e3);
    fprintf('Ground to 1st d band average splitting = %0.3f Er = %0.1f KHz \n', d_band_gaps(1), d_band_gaps(1) * Er / (hbar*2*pi) / 1e3);
    fprintf('Ground to 2nd d band average splitting = %0.3f Er = %0.1f KHz \n', d_band_gaps(2), d_band_gaps(2) * Er /(hbar*2*pi) / 1e3);
    fprintf('Ground to 3rd d band average splitting = %0.3f Er = %0.1f KHz \n', d_band_gaps(3), d_band_gaps(3) * Er / (hbar*2*pi) / 1e3);
    fprintf('Ground to next bands max gap = %0.3f Er = %0.1f KHz \n', NextHigherBandsMin, NextHigherBandsMin * Er / (hbar*2*pi) / 1e3);
end

if show_plots
    plot_band_structure(EBand, WannierFn);
end

end