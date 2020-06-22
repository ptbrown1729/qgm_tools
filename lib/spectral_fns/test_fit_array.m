xs = -4 : 4;
ys = -3 : 3;
[xx, yy] = meshgrid(xs, ys);
rs = sqrt( xx.^2 + yy.^2);
omegas = linspace(0, 10, 100);

array = zeros( length(ys), length(xs), length(omegas));
real_params = zeros( length(ys), length(xs), 4);

for ii = 1 : length(ys)
    for jj = 1 :length(xs)
        ps = [2 + rs(ii, jj), 1.111111, 1.222222, 0.113344];
        real_params(ii, jj, :) = ps;
        array(ii, jj, :) = lorentzian1D(ps, omegas) + 0.3 * rand(1, length(omegas) );
    end
end

% fitting
init_params = [14, 2, 0.5, -0.1];
fixed_params = [0, 0, 0, 0];
lbs = [min(omegas), 0, 0, -inf];
ubs = [max(omegas), inf, inf, inf];
[fit_params, std_errs, fit_ref] = fit_array_2lor(array, '', omegas, xs, ys,...
                                                 init_params, fixed_params, lbs, ubs);
                                             
% fits
omegas_interp = linspace( min(omegas), max(omegas), 300);
array_interp = zeros( length(ys), length(xs), length(omegas_interp) );
for ii = 1 : length(ys)
    for jj = 1 : length(xs)
        array_interp(ii, jj, :) = lorentzian1D(fit_params(ii, jj, :), omegas_interp);
    end
end
        
% plot
figh = plot_3d_array(array, '', xs, ys, omegas);
plot_3d_array(array_interp, '', xs, ys, omegas_interp, '', figh, 'r');

%%
dist_fn = @(x1, y1, x2, y2) min( mod( x1 - x2, 3), mod(x2 - x1, 3) );
% li = linspace(0, 3, 20);
order = dist_fn(obj.linear_index, 0, 0, 0);
init_params = [mean(obj.frqs_offsets_tunits), 1, max(obj.nr_gxmg(1, :)), 0];
fixed_params = [0, 0, 0, 1];
lbs = [min(obj.frqs_offsets_tunits), 0.5, 0, -inf];
ubs = [max(obj.frqs_offsets_tunits), inf, inf, inf];
[fp, se, refs] = fit_array_2lor(obj.nr_gxmg, '', obj.frqs_offsets_tunits, obj.linear_index, '',...
                             init_params, fixed_params, lbs, ubs, 1, dist_fn, order);

frqs_interp = linspace( min(obj.frqs_offsets_tunits), max(obj.frqs_offsets_tunits), 300);
array_interp = zeros(1, size(obj.nr_gymg, 1), length(frqs_interp) );
for ii = 1 : size(array_interp, 2)
    array_interp(1, ii, :) = lorentzian1D(fp(1, ii, :), frqs_interp);
end

% array = permute(obj.nr_gxmg, [3, 1, 2]);

figh = plot_3d_array(obj.nr_gxmg, obj.nr_gxmg_unc, obj.linear_index, '', obj.frqs_offsets_tunits, '', '', '', 'bo');
plot_3d_array(array_interp, '', obj.linear_index, '', frqs_interp, '', figh, '', 'r');
