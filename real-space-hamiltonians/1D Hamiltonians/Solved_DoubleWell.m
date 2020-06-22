[H,V,xpts]=HdoubleWell();
Neigs = 15;
[psivects,Es]=getEigs1D(H,Neigs);

Nplots = Neigs+2;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots)/Ncols);

figure('name','Double Well Eigenstates')
for ii = 1:Neigs
    subplot(Nrows,Ncols,ii)
    plot(xpts,psivects(:,ii))
end

subplot(Nrows,Ncols,Nplots-1)
plot(xpts,V(xpts))