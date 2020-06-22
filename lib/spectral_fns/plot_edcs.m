function [figh, h] = plot_edcs(mat_2d, energies, offset, figh)
% plot edcs along a 2D cut offset form each other by a certain amount.

if ~exist('offset', 'var') || isempty(offset)
    offset = 0.2 * max( mat_2d(:) );
end

offsets = offset * (1:size(mat_2d, 1));

if ~exist('figh', 'var') || isempty(figh) || ~isvalid(figh)
    figh = figure;
else
    figure(figh);
end

line_width = 2;
% plot each curve
for ii = 1:size(mat_2d, 1)

    if ii == 1
        h = plot(energies, mat_2d(ii, :) + offsets(ii), '-', 'LineWidth', line_width);
    else
        h = plot(energies, mat_2d(ii, :) + offsets(ii), '-', 'LineWidth', line_width, 'Color', h.Color);
    end
    hold on;
end


end