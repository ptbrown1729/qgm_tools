function [fig_handle] = display_high_symm_cuts(array_gxmg, array_unc_gxmg,...
                kx_gxmg, ky_gxmg, linear_index, display_name, fig_handle)
%   [fig_handle] = display_high_symm_cuts(array_gxmg, array_unc_gxmg,...
%                           kx_gxmg, ky_gxmg, linear_index, display_name, fig_handle)
% display high symmetry cuts
% If prefer gymg, can take kx_gxmg -> - ky_gymg and ky_gxmg -> kx_gxmg


if ~exist('fig_handle', 'var') || isempty(fig_handle) || ~isvalid(fig_handle)
    fig_handle = figure();
end

if isempty(array_unc_gxmg)
    array_unc_gxmg = zeros(size(array_gxmg));
end

if ~exist('display_name', 'var') || isempty(display_name)
    display_name = '';
end

[array_parts_gxmg, kx_parts_gxmg, ky_parts_gxmg] = ...
    get_high_symm_with_position( array_gxmg, kx_gxmg, ky_gxmg, linear_index);
[array_unc_parts_gxmg, ~, ~] = ...
    get_high_symm_with_position( array_unc_gxmg, kx_gxmg, ky_gxmg, linear_index);

figure(fig_handle);
nrows = 2;
ncols = 2;

% GX and GY
subplot(nrows, ncols, 1)

hold on;
ph_gx = errorbar( abs(kx_parts_gxmg{1}), array_parts_gxmg{1}, array_unc_parts_gxmg{1},...
        'DisplayName', sprintf('%s', display_name)); 
xlabel('|kx|');
title('GX');
legend('off');
legend('show');

ax = gca;
ax.YLim(1) = 0;
ax.YLim(2) = max([ax.YLim(2), 1.2 * max(array_parts_gxmg{1}) ]);
lims = ax.YLim;

% XM and YM
subplot(nrows, ncols, 2)

hold on;
ph_xm = errorbar( abs(ky_parts_gxmg{2}), array_parts_gxmg{2}, array_unc_parts_gxmg{2},...
        'DisplayName', sprintf('%s', display_name), 'Color', ph_gx.Color); 
ylim(lims);
xlabel('|ky|');
title('XM');
legend('off');
legend('show');

% MG
subplot(nrows, ncols, 3)

hold on;
ph_mg = errorbar( abs(kx_parts_gxmg{3}), array_parts_gxmg{3}, array_unc_parts_gxmg{3},...
        'DisplayName', sprintf('%s', display_name), 'Color', ph_gx.Color); 
ylim(lims);
xlabel('|kx| = |ky|');
title('MG');
legend('off');
legend('show');

subplot(nrows, ncols, 4)

hold on;
% GX
errorbar( abs(kx_parts_gxmg{1}), array_parts_gxmg{1}, array_unc_parts_gxmg{1},...
        'DisplayName', sprintf('%s GX', display_name), 'Color', ph_gx.Color); 
% MG
errorbar( sqrt( kx_parts_gxmg{3}.^2 + ky_parts_gxmg{3}.^2), array_parts_gxmg{3}, array_unc_parts_gxmg{3},...
        '--', 'DisplayName', sprintf('%s MG', display_name), 'Color', ph_gx.Color); 
xlabel('|k|');
ylim(lims);
title('GX vs. MG');
legend('off');
legend('show');
                
end