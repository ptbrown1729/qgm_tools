bin_edges_azavg = sqrt(linspace(0, 25^2, 25));
df = DataFolder({2018 4 24 57 1 1}, [], bin_edges_azavg);

bin_size = 5;
img_size = 100;
dist_grid = kron(reshape(1:(img_size/bin_size)^2, [img_size/bin_size, img_size/bin_size]), ones(bin_size));
bin_edges = -0.5:1:max(dist_grid(:))+0.5;
bin_ends = bin_edges(2:end);

[n, nsdm, r, rsdm ,npts, nI, nIsdm, nA, nAsdm, nIsdmA, nAsdmI] = ...
                df.getDensMoments(df.Occs_Stack, dist_grid, bin_edges);
% reduce to only bins with atoms in therm
bin_ends(n==0) = [];
bin_edges = [bin_edges(1), bin_ends];

% now redo
[n, nsdm, r, rsdm ,npts, nI, nIsdm, nA, nAsdm, nIsdmA, nAsdmI] = ...
    df.getDensMoments(df.Occs_Stack, dist_grid, bin_edges);
[~, ~, dist_avg, ~, ~, ~, ~] = azAvg_General(df.DistGrid, [], dist_grid, bin_edges);

            
 [nn, nnsdm, ni, nisdm, nj, njsdm, npts] = ...
    df.getCorrMoments(df.Occs_Stack, dist_grid, bin_edges, 1, 1);

[nnc, nncunc] = df.getCorrWithErr(nn, ni, nj, npts);

figure;
subplot(1, 2, 1)
errorbar(dist_avg, nnc(:, 2, 3), nncunc(:, 2, 3), 'o');
hold on;
errorbar(df.RadialPos, squeeze(df.Density_Corr_AzAvg(4, 5, :)),....
                       squeeze(df.Density_Corr_AzAvgUnc(4, 5, :)), 'rx');
grid on;
xlabel('Position (sites)');
ylabel('nn correlator');

subplot(1, 2, 2);
errorbar(n, nnc(:, 2, 3), nncunc(:, 2, 3), nncunc(:, 2, 3), nsdm, nsdm, 'o');
hold on;
errorbar(df.Occs_AzAvg, squeeze(df.Density_Corr_AzAvg(4, 5, :)),....
                       squeeze(df.Density_Corr_AzAvgUnc(4, 5, :)),...
                       squeeze(df.Density_Corr_AzAvgUnc(4, 5, :)),...
                       df.Occs_AzAvgUnc, df.Occs_AzAvgUnc, 'rx');
grid on;
xlim([-0.05, 1]);
xlabel('<n>');
ylabel('nn correlator');