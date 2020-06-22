function [t, tx, ty, tdiag, u, uoverT] = latt_2d_interpolated_py(latt_depth, efld_retro_attenuation, latt_data_file)
% [t, tx, ty, tdiag, u, uoverT] = Lattice_2D_Asymmetric_Interpolated(Depth,Attenuation)
%
% Get U/t and related lattice parameters from depth and attenuation factor.
% Uses a precomputed table of values to interpolate data because direct
% calculation is extremely slow.
%
% latt_depth:
%
% efld_retro_attenuation:
%
% latt_data_file: path to precomputed lattice data

if ~exist('latt_data_file', 'var') || isempty(latt_data_file)
    latt_data_file = 'lattice_data.txt';
end

%Read data.
latt_data = dlmread(latt_data_file, '\t' , 1, 0);

%
depths = latt_data(:, 1);
efield_attn = latt_data(:, 2);

%Reshape these into meshgrid like shapes. This makes interpolation much
%faster.
n_depths = length(unique(depths));

% if scan retro then depths
% depths_grid = reshape(depths, [length(depths) / n_depths, n_depths]);
% retro_atten_grid = reshape(efield_attn, [length(depths) / n_depths, n_depths]);

depths_grid = reshape(depths, [n_depths, length(depths) / n_depths]);
retro_atten_grid = reshape(efield_attn, [n_depths, length(depths) / n_depths]);

%reshape the data to have a similar shape, where the third axis is now the
%various different lattice parameters.
latt_data_grid = reshape(latt_data, [n_depths, length(depths) / n_depths, size(latt_data, 2)]);

txs = latt_data_grid(:, :, 5);
tys = latt_data_grid(:, :, 6);
tds = latt_data_grid(:, :, 7);
us = latt_data_grid(:, :, 8);

%interpolating functions
tx_fn = @(depth, r) interp2(retro_atten_grid, depths_grid, txs, r, depth);
ty_fn = @(depth, r) interp2(retro_atten_grid, depths_grid, tys, r, depth);
td_fn = @(depth, r) interp2(retro_atten_grid, depths_grid, tds, r, depth);
u_fn = @(depth, r) interp2(retro_atten_grid, depths_grid, us, r, depth);

tx = tx_fn(latt_depth, efld_retro_attenuation);
ty = ty_fn(latt_depth, efld_retro_attenuation);
tdiag = td_fn(latt_depth, efld_retro_attenuation);
u = u_fn(latt_depth, efld_retro_attenuation);
uoverT = u./t;

end