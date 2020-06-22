[H,V,xpts,omegaHarm]=Hgaussian1D(1,120e-6);
Neigs = 50;
[psivects,Es]=getEigs1D(H,Neigs);

Nplots = Neigs+2;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots)/Ncols);

figure('name','Gaussian1D')
for ii = 1:Neigs
    subplot(Nrows,Ncols,ii)
    plot(xpts,psivects(:,ii))
end

subplot(Nrows,Ncols,Nplots-1)
plot(xpts,V(xpts))

subplot(Nrows,Ncols,Nplots)
plot(repmat(Es,2)'/hbar/omegaHarm-0.5)