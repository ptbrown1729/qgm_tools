%for testing 1D SHO Hamiltonian compared with anlytic result.

%run('..\constants.m')
%constants
mLi = 9.9883e-27;
hbar = 1.0546e-34;
c = 299792458;
abohr = 5.2918e-11;

Neigs = 15;

omega=2*pi*3e3;
[H, V , xpts] = Hsho1D(omega);
[psivects, Es] = getEigs1D(H, Neigs);

Nplots = Neigs+1;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots) / Ncols);

%plot eigenstates, and compare to analytic result
figure('name','1D SHO, Analytic Vs. Numeric')
for ii = 1 : (Nplots - 1)
    
subplot(Nrows,Ncols,ii)
plot(xpts, psivects(:,ii) / max(abs(psivects(:,ii))) / sign(psivects(1,ii)), 'b.')
hold on;
plot(xpts, getPsi1DSHO(xpts, ii-1) / max(abs(getPsi1DSHO(xpts, ii-1))) / sign(max(getPsi1DSHO(xpts(1), ii-1))))
end
subplot(Nrows,Ncols,Nplots)
plot(xpts,V(xpts))
hold off;

%ground state, fourth power expectation value.
Gstate_4th_analytic = sqrt(mLi*omega/pi/hbar/2);
dx = xpts(2) - xpts(1);
Gstate_4th_numeric = sum((psivects(:,1).*conj(psivects(:,1))).^2)*dx/dx^2;