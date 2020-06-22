function [wanns,xpositions] = Wannier1D(H,xpts,Nlowband)
%Finds Wannier functions for 1D case (with no band degeneracies?) by
%diagonalizing position operator projected on to lowest band.
%This method does not generalize to higher dimensions (?)
%Based on Erich Mueller presentation for Topological and Strongly Correlated
%Phases in Cold Atoms conference at Princeton 04/29/15.

[psivects,~]=getEigs1D(H,Nlowband);
Nx = length(H);
%build projector on to lowest band
ProjForm = sparse(1:Nx,ones(1,Nx),psivects(:,1),Nx,Nx);
for ii=2:Nlowband
    ProjForm = ProjForm + sparse(1:Nx,ii*ones(1,Nx),psivects(:,ii),Nx,Nx);
end
%ProjForm(:,1:Nlowband) = psivects(:,1:Nlowband);
P = ProjForm*ProjForm';
X = sparse(1:Nx,1:Nx,xpts);
Xlowband = P*X*P;
[wanns,xpositions] = eigs(Xlowband,Nlowband);
wanns=wanns(:,1:Nlowband);
[xpositions,indices] = sort(diag(xpositions));
wanns = wanns(:,indices);

%fix sign so function is positive.
[~,maxIndices] = max(abs(wanns));
columnIndices = 1:size(wanns,2);
signs = sign(wanns(sub2ind(size(wanns),maxIndices,columnIndices)));
wanns = wanns.*kron(ones(size(wanns,1),1),signs);
end

