function [psivects,Es] = getEigs1D(H,Neigs)
%[psivects,Es] = getEigs1D(H,Neigs)
%Finds the Neigs smalled magnitude eigenvalues of H
%Note that the potential V should be everywhere greater than zero for this
%to work...since otherwise the ground state energy becomes negative, and
%may no longer have the smallest magnitude.

%get eigenvalues and eigenvectors
[psivects,Emat] = eigs(H,Neigs,'sm');
[Es,position] = sort(diag(Emat));
psivects=psivects(:,position);
end

