%E(k) = 2*(1 - cos(k))
dos = @(e) 1/(2*pi) * 1 ./ sqrt(1 - (1 - e/2).^2);
cumulative_states = @(e) 1/(pi) * ( asin(1) - asin(1 - e/2));
f = @(e, T) 1./ (exp(e./T) + 1);

%function of (mu,T)
n = @(mu, T)integral(@(e) dos(e) .* f(e - mu, T), 0, 4, 'ArrayValued',true);
chi = @(mu, T) integral(@(e) dos(e) .* f(e - mu, T).^2 .* exp( (e - mu)./ T) ./ T, 0, 4, 'ArrayValued',true);

%convert to functions of (n,T)
%desired values
n_gridpts = 300;
mus_all = linspace(-10, 10, n_gridpts);
Ts_all = linspace(0, 10, n_gridpts);
ns_all = linspace(0.02, 1, n_gridpts);
%grids
[mu_grid, T_grid] = meshgrid(mus_all, Ts_all);
n_on_mu_grid = n(mu_grid, T_grid);
mu = scatteredInterpolant(n_on_mu_grid(:), T_grid(:), mu_grid(:));
[n_grid, ~] = meshgrid(ns_all, Ts_all);


%display Fermi gas data
es = linspace(0,4,100);

figure;
subplot(3,2,1);
plot(es, dos(es));
xlabel('Energy (t)');
ylabel('dos(e)');
grid on;

subplot(3,2,2);
plot(es, cumulative_states(es));
xlabel('Energy (t)');
ylabel('cumulative states (e)');
grid on;

subplot(3,2,3)
plot(es, f(es-3, 1));
xlabel('Energy (t)');
ylabel('f(e)');
grid on;

subplot(3,2,4);
imagesc(mu(n_grid, T_grid));
axis equal; axis image;
title('\mu');
xlabel('n');
ylabel('T');

subplot(3,2,5);
imagesc(chi(mu(n_grid, T_grid), T_grid));
axis equal; axis image;
title('dn/dmu');
xlabel('n');
ylabel('T');

%plot compressibility vs. T
subplot(3,2,6)
ts = logspace(-2,2,1000);
ns = 0.83 * ones(size(ts));
loglog(ts, 1./chi(mu(ns, ts), ts));
grid on;
xlabel('Temp (t)');
ylabel('dn/dmu');
title(sprintf('n = %0.2f', ns(1)));

suptitle('1D Single-Band Non-Interacting Fermi Gas');