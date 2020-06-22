omegas = linspace(0, 15, 100);
beta = 1;
mu = -2;

chis = zeros(50, 50, length(omegas));
for ii = 1:length(omegas)
    [kx, ky, chi] = fg_compress(omegas(ii), beta, mu);
    chis(:, :, ii) = fftshift(chi);
end
chis(chis < 0) = 0;
kx = fftshift(kx);
ky = fftshift(ky);

kx( abs(kx - pi) < 1e-15) = pi;
ky( abs(ky - pi) < 1e-15) = pi;
kx( kx >= pi ) = kx( kx >= pi) - 2*pi;
ky( ky >= pi ) = ky( ky >= pi) - 2*pi;

figure;
index1 = 26;
index2 = 27;

plot(omegas, squeeze(chis(index1, index2, :)));
xlabel('omegas (t)');
ylabel('chi(k, omega) (1/ta^2)');
ttl = sprintf('kx = %0.2f pi, ky = %0.2f pi', kx(index1, index2)/pi, ky(index1, index2)/pi);
title(ttl);
xlim([min(omegas), max(omegas)]);

%%
T = linspace(0.1,10,100);
mu = 0;

chi0 = zeros(1, length(T));
chi1 = zeros(1, length(T));
for ii = 1:length(chi0)
    chi0(ii) = 2*fg_compress_static(1/T(ii), mu);
    chi1(ii) = FG.chifn_mu_dmu_T(mu,0,T(ii));
end

figure;
plot(T, chi0, 'b.');
hold on;
plot(T, chi1, 'r.');
xlabel('T (t)');
ylabel('\chi_c (1/ta^2)');
grid on;


%% test


delta_mu = 0.01;
temp = 0.25;
diff_delta_mu = (fg_density(1/temp, delta_mu) - fg_density(1/temp, -delta_mu)) / (2 * delta_mu)

fg_compress_static(1/temp, 0)
FG.chifn_mu_dmu_T(0, 0, temp)/2


%%
mu = 0;
T = linspace(0.3, 15, 200);
nsites = 50:20:200;

chi = zeros(length(T), length(nsites));

for ii = 1:length(T)
    for jj = 1:length(nsites)
        chi(ii, jj) = fg_compress_static(1/T(ii), mu, nsites(jj));
    end
end

figure;
loglog(T, chi, '.');
xlabel('T(t)');
ylabel('\chi 1/(ta^2)');
grid on;

leg = {};
for ii = 1:length(nsites)
    leg{ii} = sprintf('nsites = %d', nsites(ii));
end
legend(leg);

