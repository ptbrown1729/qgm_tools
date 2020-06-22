%test real space and k-space diagonalizations of lattice.
run('constants.m')

%lattice parameters
lambda=1064e-9;
k=2*pi/lambda;
Er=hbar^2*k^2/(2*mLi);

V = @(x) Er*cos(2*k*x);

%Get x-points.
Nbrill = 201;
%initialize grid
%dx should be a fraction of a.
a = lambda/2;
dx = a/50; 
xstart = -Nbrill*a/2;
xend = Nbrill*a/2-dx;
% xstart = 0;
% xend = Nbrill*a-dx;
Nx = round((xend-xstart)/dx+1);
xpts = xstart:dx:xend;


H=Hlattice(V,lambda/2);
[psivects,EsReal]=getDiagH(H,201);

[qvects,EsKspace,psiKs,recpvects] = OneDBand(V,lambda/2);

%Test one: energies
subplot(2,2,1)
plotxs = zeros(201,1);
plot(plotxs,EsReal/Er,'r.')
hold on;
plot(plotxs+0.1,EsKspace(:,1)/Er,'b.')
xlim([-0.1,0.2])

%Test two: check that get same wavefunctions
%q=0 ground band wavefunction from k-space diagonalization
subplot(2,2,2)
qvZeros = psiKs(:,:,101);
qvZerosGroundBand = qvZeros(:,1);
qvZerosGroundBandReal = getPsiRealFromPsiK(xpts,qvZerosGroundBand,recpvects,qvects(101));
qvZerosGroundBandReal = qvZerosGroundBandReal/sqrt(qvZerosGroundBandReal*qvZerosGroundBandReal');
plot(xpts,qvZerosGroundBandReal,'r')
hold on
plot(xpts,psivects(:,1),'b.')

%next q value. Since degenerate values, this will not necessarily be a
%Bloch wavefunction.
psiRone = psivects(:,2);
psiRtwo = psivects(:,3);
psiKone = psiKs(:,1,100);
psiKoneR = getPsiRealFromPsiK(xpts,psiKone,recpvects,qvects(100));
psiKoneR = psiKoneR/sqrt(psiKoneR*psiKoneR');
psiKtwo = psiKs(:,1,102);
psiKtwoR = getPsiRealFromPsiK(xpts,psiKtwo,recpvects,qvects(102));
psiKtwoR = psiKtwoR/sqrt(psiKtwoR*psiKtwoR');
test1 = abs(psiKoneR*psiRone)^2+abs(psiKtwoR*psiRone)^2
test2 = abs(psiKoneR*psiRtwo)^2+abs(psiKtwoR*psiRtwo)^2


%Test 3: calculate Wannier states using both methods.
% subplot(2,2,3)
% [WannsReal,xpositions] = Wannier1D(H,xpts,51);
% hold all
% for ii=Nbrill+3:Nbrill-3
%     plot(xpts-xpositions(ii),WannsReal(:,ii))
% end

subplot(2,2,4)
Nband = 1;
[wannFromK,xpositionFromK] = getWannier1DFromBloch(xpts,psiKs,recpvects,qvects,Nband);

plot(xpts-xpositionFromK,wannFromK)






% for ii=1:10
% subplot(2,5,ii)
% plot(0:lambda/100:lambda*51/2-lambda/100,psivects(:,ii)/max(abs(psivects(:,ii)))/sign(psivects(1,ii)))
% end
% hold off
% figure
% plot(0.1*ones(1,51),Es/Er,'r.')