function [psiT] = evolveTimeIndH(startpsi,eigenstates,energies,t)
%[psiT] = evolveTimeIndH(startpsi,eigenstates,energies,t
%psis are column vectors. eigenstates is a length x number
%matrix of eigenstates

%conj(transpose(BM))*(inital basis) gives element in new basis.
BasisChangeMat = conj(transpose(eigenstates));
startpsiEBasis = BasisChangeMat*startpsi;
hbar = 1.0546e-34;
psiTEBasis =  startpsiEBasis.*exp(-1i*energies/hbar*t);
psiT = eigenstates*psiTEBasis;

end