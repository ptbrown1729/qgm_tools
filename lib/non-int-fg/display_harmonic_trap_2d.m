function fig_handle = display_harmonic_trap_2d(mu_struct, axes)
% display results from fit_harmonic_trap_2d function


if ~exist('axes', 'var')
    fig_handle = figure;
    axes  = gca;
else
    fig_handle = axes.Parent;
end

if is_meshgrid(mu_struct.xpos, mu_struct.ypos)
% if we have a meshgrid, display using imagesc 
    mu_and_diff_gauss = horzcat(mu_struct.mus, mu_struct.fn_mu_gauss(mu_struct.xpos, mu_struct.ypos) - mu_struct.mus);
    mu_and_diff_parab = horzcat(mu_struct.mus, mu_struct.fn_mu_parab(mu_struct.xpos, mu_struct.ypos) - mu_struct.mus);
    pic = vertcat(mu_and_diff_gauss, mu_and_diff_parab);
    
%     subplot(1, 2, 1);
    imagesc(axes, pic);
    colorbar;
    axes.XTickLabel = '';
    axes.YTickLabel = '';
    axis equal;
    axis image;
    
%     clims = caxis;
%     cmax = clims(2);
%     
%     subplot(1, 2, 2)
%     imagesc( , [-cmax, cmax]);
%     colorbar;
%     axis equal;
%     axis image;
    
    color_map = get_color_map(100, [0, 0, 1], [1, 1, 1], [1, 0, 0]);
    colormap(color_map);
    
else
    % otherwise display as 3D plot

    [xx_interp, yy_interp] = meshgrid( min(mu_struct.ypos(:)) : max(mu_struct.ypos(:)),...
                                       min(mu_struct.xpos(:)) : max(mu_struct.xpos(:)));
    
    mesh(axes, xx_interp, yy_interp,...
         mu_struct.fn_mu_gauss(xx_interp, yy_interp));
    hold on;
    scatter3(axes, mu_struct.xpos, mu_struct.ypos, mu_struct.mus,...
             300, mu_struct.mus, '.' );
end

% print information about fits
gauss_str = sprintf('Wx = %0.2f(%.0f), Wy = %0.2f(%.0f) \n wx = (2pi) %0.2f(%.0f), wy = (2pi) %.2f(%.0f), theta = %0.0f(%.0f) deg',...
        mu_struct.fit_params_mu_gauss(3),...
        mu_struct.fit_err_mu_gauss(3),...
        mu_struct.fit_params_mu_gauss(4),...
        mu_struct.fit_err_mu_gauss(4),...
        mu_struct.omegax_gauss_sqrt_toverh / (2*pi),...
        mu_struct.omegax_gauss_sqrt_toverh_unc / (2*pi) * 1e2,...
        mu_struct.omegay_gauss_sqrt_toverh / (2*pi),...
        mu_struct.omegay_gauss_sqrt_toverh_unc / (2*pi) * 1e2,...
        mu_struct.fit_params_mu_gauss(6) * 180 / pi,...
        mu_struct.fit_err_mu_gauss(6)  * 180 / pi );
    
parab_str = sprintf('wx parab = (2pi)%0.2f (%.0f), wy parab = (2pi)%0.2f(%.0f), theta = %.0f(%.0f)deg',...
        mu_struct.omegax_parab_sqrt_toverh / (2*pi),...
        mu_struct.omegax_parab_sqrt_toverh_unc / (2*pi) * 1e2,...
        mu_struct.omegay_parab_sqrt_toverh / (2*pi),...
        mu_struct.omegay_parab_sqrt_toverh_unc / (2*pi) * 1e2,...
        mu_struct.fit_params_mu_parab(5) * 180 / pi,...
        mu_struct.fit_err_mu_parab(5) * 180 / pi);
    
ttl = sprintf('waist (a), w sqrt(t/h), %s\n%s', gauss_str, parab_str);
title(ttl);


end

function is_grid = is_meshgrid(xx, yy)

is_grid = 0;

if ~ismatrix(xx) || ~ismatrix(yy)
    return;
end

if size(xx) ~= size(yy)
    return;
end

for jj = 1:size(xx, 1)
    if any(xx(:, jj) ~= xx(1, jj))
        return;
    end
end

for ii = 1 : size(yy, 1)
    if any( yy(ii, :) ~= yy(ii, 1) )
        return;
    end
end

is_grid = 1;

end