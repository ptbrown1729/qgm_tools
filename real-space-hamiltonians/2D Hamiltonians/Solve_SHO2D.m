Neigs = 15;
[H,V,xpts,ypts] = Hsho2D();
[psivects,Es] = getEigs2D(H,Neigs);

Nplots = Neigs+2;
Ncols = ceil(sqrt(Nplots));
Nrows = ceil((Nplots)/Ncols);

figure('name','2D SHO')
eigvectorStack = zeros(length(ypts),length(xpts),Neigs);
for ii = 1:Nplots-2
    eigvectorStack(:,:,ii) = reshape(psivects(:,ii),[length(ypts),length(xpts)]);
    subplot(Nrows,Ncols,ii)
    imagesc(eigvectorStack(:,:,ii));
    axis equal;
    axis image;

end

subplot(Nrows,Ncols,Nplots-1)
[X,Y] = meshgrid(xpts,ypts);
imagesc(V(X,Y));
axis equal;
axis image;