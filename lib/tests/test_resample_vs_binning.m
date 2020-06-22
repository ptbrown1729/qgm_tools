% Test for difference in binning versus resampling functions.

n = 200;
nbin = 4;
img = rand(n, n);
unc = rand(n, n);

[b, bu, bsd, bsem] = binImg(img, nbin, nbin, 'Normal', unc);
[r, dx, dky, rsd, rsem, ru] = get_resampled_img(img, [n/nbin, n/nbin], '', '', 'mean', unc);

max_diff = max(max(abs(b -r)));
max_diff_sd = max(max(abs(bsd - rsd)));
max_diff_sem = max(max(abs(bsem - rsem)));
max_diff_unc = max(max(abs(bu - ru)));

fprintf('testing binning versus resampling\n');
fprintf('max differince in value is %0.3e\n', max_diff);
fprintf('max differince in std is %0.3e\n', max_diff_sd);
fprintf('max differince in sem is %0.3e\n', max_diff_sem);
fprintf('max difference in unc is %0.3e\n', max_diff_unc);
fprintf('\n');