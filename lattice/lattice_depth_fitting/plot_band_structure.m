function [fig_handle] = plot_band_structure(energy_bands, wannier_fn, n_bands_to_plot)
% Display band structure and wannier function produced by
% Lattice_2D_Asymmetric
%
% EBand: array of energy bands, with size qy_pts x qx_pts x n_bands
%
% WannierFn: wannier function of the ground band

if ~exist('wannier_fn', 'var') || isempty(wannier_fn)
    wannier_fn = 0;
end

if ~exist('n_bands_to_plot', 'var') || isempty(n_bands_to_plot)
    n_bands_to_plot = 6;
end

hbar = 1.055e-34; 
m = 6.015*1.66e-27;
a = 1064e-9 / (2^0.5); 
e_recoil = (pi * hbar)^2 / (2 * m * a^2);

% get bands along Gamma-M-X-Gamma
e_bands_quad = get_quadrant(permute(energy_bands(1 : n_bands_to_plot, :, :), [2, 3, 1]), 'northeast');
[ebands_gxmg, ~, ~, ~, linear_index] = get_high_symm_cut(e_bands_quad, 'lowerleft');

%plot information about lowest bands
% figtitle = sprintf('Latt Depth = %0.1f Er', s);
fig_handle = figure();

% plot bands along GMXG path
subplot(1, 3, 1)
e_gxmg_offset = ebands_gxmg - min(energy_bands(:));
plot(linear_index, e_gxmg_offset, '.-');
title('Band structure, high symmetry pts')
ylabel('Energy (Er)')
xticks([0, 1, 2, 3]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'}); 

% 3D band plot
subplot(1, 3, 2)
for ii = 1 : n_bands_to_plot
    surf(squeeze(energy_bands(ii,:,:) - min(energy_bands(:))), 'EdgeColor', 'none');
    hold on;
end
title('2D Band Structure')
zlabel('Energy (Er)')

% plot wannier functions
subplot(2, 3, 3);
imagesc(real(wannier_fn))
axis equal;
axis image;
title('Wannier Function, Real Part')

subplot(2, 3, 6);
imagesc(imag(wannier_fn))
axis equal;
axis image;
title('Wannier Imaginary, Real Part') 

end