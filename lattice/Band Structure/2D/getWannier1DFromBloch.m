function [wanns,R] = getWannier1DFromBloch(xpts,psiKs,recpvects,qvects,Nband)
%psiKs is a Kvects x Kvects x Nbrill matrix. psiK x band index x qvalue
%Nbrill = length(qvects)
%Nband is the band two get Wannier functions from.

%permute to get rid of extra dimension.
psiKs = permute(psiKs(:,Nband,:),[1,3,2]);

psiKsReal = zeros(length(xpts),length(qvects));
for ii=1:length(qvects)
psiKsReal(:,ii) = getPsiRealFromPsiK(xpts,psiKs(:,ii),recpvects,qvects(ii));
end

Rpts = 2*pi./recpvects;
Rpts(isinf(Rpts))=0;
R = Rpts(floor(length(Rpts)/2));

%wann_R(x) = sum_q exp(-1i*q*R) psi_q(r)
%create 3D array, then sum over one of the dimensions.
% [qvals,~,Rvals] = meshgrid(qvects,xpts,Rpts);
% tempWanns = exp(-1i*Rvals.*qvals).*repmat(psiKsReal,[1,1,length(Rpts)]);
% wanns = permute(sum(tempWanns,2),[1,3,2]);

[qvals,~] = meshgrid(qvects,xpts);
tempwanns = exp(-1i*R*qvals).*psiKsReal;
wanns = sum(tempwanns,2);
