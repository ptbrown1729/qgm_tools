[H,V,xpts]=HdoubleWell();
Neigs = 100;
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

h=6.626e-34;
tf = h/Es(5)*20.56353;
Nts = 50;
ts = linspace(0,tf,Nts);

Nplots = Nts;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots)/Ncols); 


InitialState = circshift(psivects(:,1),[200,0]);
psiTs = zeros(length(InitialState),Nts);
figure('name','Time Evolution')
for ii = 1:Nts
    psiTs(:,ii) = evolveTimeIndH(InitialState,psivects,Es,ts(ii));
    subplot(Nrows,Ncols,ii);
    plot(abs(psiTs(:,ii)).^2);
end
psiTsMags = abs(psiTs).^2;


