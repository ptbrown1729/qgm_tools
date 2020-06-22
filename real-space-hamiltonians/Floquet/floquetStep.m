function [eigvects,eigvals]=floquetStep(H,T)
%Floquet calculation using incremental evolution of operator
%H is Hamiltonian matrix as a function of time given in some basis.
%H(t+T)=H(t). H is T-periodic.

run('constants.m');

steps=100;
dt=T/steps;
Hexp=@(t) expm(-1i*H(t)*dt/hbar);

dim = length(H(0));
Uinit=sparse(1:dim,1:dim,ones(1,dim));
%Uinit = eye(size(H(0));
U=Uinit;
for ii=1:steps
    U=Hexp(ii*dt)*U;
end

HF=1i*hbar*logm(U)/T;

if issparse(HF)==1
    [eigvects,Es]=eigs(HF,5,'sm');
    eigvals=real(diag(Es));
else
    [eigvects,Es]=eig(HF);
    EsFl=mod(real(diag(Es)),hbar*2*pi/T);
    [eigvals,sortedIndices]=sort(EsFl);
    eigvects = eigvects(:,sortedIndices);
end
end
