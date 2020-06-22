function [psivects,Es] = getEigs2D(H,Neigs)
%[psivects,Es] = getEigs2D(H,Neigs)
%psivects is an N_H by Neigs vector, i.e. each column is a single eigenvector

%get eigenvalues and eigenvectors
[psivects,Emat] = eigs(H,Neigs,'sm');
[Es,position] = sort(diag(Emat));
psivects=psivects(:,position);
end

