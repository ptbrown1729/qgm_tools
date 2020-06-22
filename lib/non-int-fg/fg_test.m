mu = 0.3;
T = linspace(0.1, 1, 200);
nsites = 10:20:200;

% density test
% susceptibility test

n = zeros(length(T), length(nsites));
chi = zeros(length(T), length(nsites));

for ii = 1:length(T)
    for jj = 1:length(nsites)
        n(ii, jj) = fg_density(1/T(ii), mu, nsites(jj));
        chi(ii, jj) = fg_compress_static(1/T(ii), mu, nsites(jj));
    end
end

leg = {};
for ii = 1:length(nsites)
    leg{ii} = sprintf('nsites = %d', nsites(ii));
end

% density
fig1 = figure;
loglog(T, n, '.');
xlabel('T(t)');
ylabel('n 1/a^2');
grid on;
legend(leg);

% susceptibility
fig2 = figure;
loglog(T, chi, '.');
xlabel('T(t)');
ylabel('\chi 1/(ta^2)');
grid on;
legend(leg);

%% test mu solver
T = linspace(0.1, 1, 100);
n = 0.3;
nsites = 10:20:50;

mu = zeros(length(T), length(nsites));

for ii = 1:length(T)
    for jj = 1:length(nsites)
        mu(ii, jj) = fg_mu(1/T(ii), n, nsites(jj));
    end
end

leg = {};
for ii = 1:length(nsites)
    leg{ii} = sprintf('nsites = %d', nsites(ii));
end

% density
fig1 = figure;
loglog(T, mu, '.');
xlabel('T(t)');
ylabel('mu (t)');
grid on;
legend(leg);
