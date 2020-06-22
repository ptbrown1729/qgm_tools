Neigs = 1;

%Fundamental constants
hbar = 1.0546e-34;
h = hbar*(2*pi);
mLi = 9.9883e-27;
abohr = 5.2918e-11;

%trapping frequencies, scattering lengths...
a = 2700*abohr;
OmegaRad = 2*pi*75e3;
OmegaAxial = 2*pi*20e3;
acharRad = sqrt(hbar/(mLi*OmegaRad));
acharAxial = sqrt(hbar/(mLi*2*pi*OmegaAxial));

%solve interacting problem 1D
[H,V,x1pts,x2pts] = H_TwoParticles_DeltaFunction_HarmonicTrap([a*1/(2*pi*acharRad^2),OmegaAxial]);
[psivectsInt,Es] = getEigs2D(H,Neigs);

%non interacting problem
H = H_TwoParticles_DeltaFunction_HarmonicTrap([0,OmegaAxial]);
[psivectsNonInt,~] = getEigs2D(H,Neigs);

%non-int SHO, theory
psiAnalytic = @(x) (mLi*OmegaAxial/pi/hbar)^(1/4)*exp(-mLi*OmegaAxial/2/hbar*x.^2);
psiAnalytic2D = @(x,y) psiAnalytic(x).*psiAnalytic(y);
[X1,X2] = meshgrid(x1pts,x2pts);
psiAnalyticMat = psiAnalytic2D(X1,X2); 

Nplots = 4;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots)/Ncols);

%figure('name',' SHO')
eigvectorStack = zeros(length(x2pts),length(x1pts),Neigs);
eigvectorStackNonInt = zeros(length(x2pts),length(x1pts),Neigs);
eigvectorDifference = zeros(length(x2pts),length(x1pts),Neigs);

for ii = 1:Neigs
eigvectorStack(:,:,ii) = reshape(psivectsInt(:,ii),[length(x2pts),length(x1pts)]);
eigvectorStackNonInt(:,:,ii) = reshape(psivectsNonInt(:,ii),[length(x2pts),length(x1pts)]);
eigvectorDifference(:,:,ii) = reshape(psivectsInt(:,ii)-psivectsNonInt(:,ii),[length(x2pts),length(x1pts)]);

figure('name',sprintf('Eigenstate %i',ii))
subplot(Nrows,Ncols,1)
%imagesc(eigvectorStack(:,:,ii));
imagesc(eigvectorDifference(:,:,ii));
axis equal;
axis image;
title('Difference Int/NonInt');
colorbar;

subplot(Nrows,Ncols,2)
imagesc(eigvectorStack(:,:,ii))
axis equal;
axis image;
title('Interacting Eigenstate')
ColorLims = caxis;

subplot(Nrows,Ncols,3)
imagesc(eigvectorStackNonInt(:,:,ii),ColorLims)
axis equal
axis image;
title('NonInteracting Eigenstate')

subplot(Nrows,Ncols,4)
imagesc(psiAnalyticMat,ColorLims);
axis equal;
axis image;
title('Analytic Result');


end






g = 4*pi*hbar^2*a/mLi;
dx1 = x1pts(2)-x1pts(1);
dx2 = x2pts(2)-x2pts(1);
%want something other than this...
psiNonInt_x1_x1 = diag(eigvectorStackNonInt(:,:,1));
UNonInt = sum((psiNonInt_x1_x1.*conj(psiNonInt_x1_x1)).^2)*dx1*dx2/(dx1^2*dx2^2);
%UNonInt = g*sum(sum((psivectsNonInt(:,:,1).*conj(psivectsNonInt(:,:,1))).^2))*dx1*dx2/(dx1^2*dx2^2);%/sum(sum(psivectsNonInt(:,:,1).*conj(psivectsNonInt(:,:,1)))); 
UNonIntThry = sqrt(mLi*OmegaAxial/pi/hbar/2);
UInt = sum(sum((eigvectorStack(:,:,1).*conj(eigvectorStack(:,:,1))).^2))*dx1*dx2/(dx1^2*dx2^2); 


overlap = sum(sum(psivectsNonInt(:,:,1).*conj(psivectsInt(:,:,1))))/sqrt(sum(sum(psivectsInt(:,:,1).*conj(psivectsInt(:,:,1))))*sum(sum(psivectsNonInt(:,:,1).*conj(psivectsNonInt(:,:,1)))));
fprintf('Assuming standard SHO eigenfunction, U = %0.1f \n',UNonInt/h);
fprintf('Modified eigenfunction, U = %0.1f \n',UInt/h);
fprintf('Umod/U = %0.2f \n',UInt/UNonInt);
fprintf('Overlap of this function with SHO ground state is = %0.5f \n',overlap) 
% subplot(Nrows,Ncols,Nplots-1)
% [X,Y] = meshgrid(x1pts,x2pts);
% imagesc(V(X,Y));
% axis equal;
% axis image;