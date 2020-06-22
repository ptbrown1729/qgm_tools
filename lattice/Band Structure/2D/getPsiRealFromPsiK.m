function [psiX] = getPsiRealFromPsiK(x,psiK,RecpVects,q)
%x is the point the wave function is to be evaluated at.
%psiK is the K-space wavefunction. A matrix with each entry corresponding
%to RecpVects
%q is the associated value in the first Brillouin zone.

if ~isequal(size(psiK),size(RecpVects))
    psiX=0;
else

%function is vectorized.
RecpVectsMatr=zeros(size(x,1),size(x,2),length(RecpVects));
for ii=1:length(RecpVects)
    RecpVectsMatr(:,:,ii)=psiK(ii)*exp(-1i*RecpVects(ii)*x);
end

psiX=sum(RecpVectsMatr,3).*exp(1i*q*x);
psiX = psiX/sqrt(psiX*psiX');

end
end

