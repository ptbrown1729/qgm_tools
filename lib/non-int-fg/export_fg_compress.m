fname_out = 'non_int_fg_compress.dat';

temps = linspace(0.1, 15, 200);
nsites = 400;
n = 0.825/2;

mu = zeros(length(temps), 1);
chi = zeros(length(temps), 1);
for ii = 1:length(temps)
    mu(ii) = fg_mu(1/temps(ii), n, nsites);
    chi(ii) = fg_compress_static(1/temps(ii), mu(ii), nsites);
end

% figure;
plot(temps, chi);
ylabel('\chi_c 1/(at^2)');
xlabel('T (t)');
grid on;

comment = sprintf('#n = %0.3f, nsites = %d', n, nsites);
names_cell = {'temp(t)', 'mu(t)', 'chi_c(1/ta^2)'};
dat_array = horzcat(temps', mu, chi);
save_data_file(dat_array, names_cell, '\t', comment, pwd, fname_out, 0);

