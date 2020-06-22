function [t, tx, ty, tdiag, u, uoverT] = ...
    Lattice_2D_Asymmetric_Interpolated(latt_depth, efld_retro_attenuation, latt_data_file)
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
    latt_data_file = 'LatticeBandData_Interpolation.txt';
end

%Read data.
latt_data = dlmread(latt_data_file, ',' , 1, 0);
depths = latt_data(:, 1);
efield_attn = latt_data(:, 2);

%Reshape these into meshgrid like shapes. This makes interpolation much
%faster. Note that assumption here is first scanned attenuation at each
%depth, then incremented depth
n_depths = length(unique(depths));
depths_grid = reshape(depths, [length(depths) / n_depths, n_depths]);
retro_atten_grid = reshape(efield_attn, [length(depths) / n_depths, n_depths]);

%reshape the data to have a similar shape, where the third axis is now the
%various different lattice parameters.
latt_data_grid = reshape(latt_data, [length(depths) / n_depths, n_depths, size(latt_data, 2)]);
ts = latt_data_grid(:, :, 6);
txs = latt_data_grid(:, :, 7);
tys = latt_data_grid(:, :, 8);
tds = latt_data_grid(:, :, 9);
us = latt_data_grid(:, :, 10);

%interpolating functions
t_fn = @(Depth,Atten) interp2(depths_grid, retro_atten_grid, ts, Depth, Atten);
tx_fn = @(Depth,Atten) interp2(depths_grid, retro_atten_grid, txs, Depth, Atten);
ty_fn = @(Depth,Atten) interp2(depths_grid, retro_atten_grid, tys, Depth, Atten);
td_fn = @(Depth,Atten) interp2(depths_grid, retro_atten_grid, tds, Depth, Atten);
u_fn = @(Depth,Atten) interp2(depths_grid, retro_atten_grid, us, Depth, Atten);

t = t_fn(latt_depth, efld_retro_attenuation);
tx = tx_fn(latt_depth, efld_retro_attenuation);
ty = ty_fn(latt_depth, efld_retro_attenuation);
tdiag = td_fn(latt_depth, efld_retro_attenuation);
u = u_fn(latt_depth, efld_retro_attenuation);
uoverT = u./t;

end