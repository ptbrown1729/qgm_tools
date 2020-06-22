sizex = 300;
sizey = 300;

[X,Y] = meshgrid(1:sizex, 1:sizey);

p = [150, 150, 30, 20, 1, pi/4, 0.3];
m = gaussian2D(p, X, Y) + gaussian2D([50, 50, 25, 25, 1, 0, 0] , X, Y) + 1*(rand(sizey, sizex) - 0.5);

initp = [127, 125, 8, 10, 0.5, 0, 0];
[fp, pn, ffh, se] = fit2D(X, Y, m, [], {'gaussian2D'}, initp, [], [], [], [], 'fminunc');

figure;
subplot(1,2,1)
imagesc(m); axis equal; axis image;

subplot(1,2,2)
imagesc(ffh(X,Y)); axis equal; axis image;

suptitle(sprintf(...
    'Cx = %0.2f -> %0.2f +/ %0.2f, Cy = %0.2f -> %0.2f +/- %0.2f\n Sx = %0.2f -> %0.2f +/- %0.2f, Sy = %0.2f -> %0.2f +/- %0.2f \n A = %0.2f -> %0.2f +/- %0.2f, Bg = %0.2f -> %0.2f +/- %0.2f\n Theta = %0.2f -> %0.2f +/- %0.2f',...
    p(1), fp(1), se(1), p(2), fp(2), se(2), ...
    p(3), fp(3), se(3), p(4), fp(4), se(4), ...
    p(5), fp(5), se(5), p(7), fp(7), se(7), ...
    p(6), fp(6), se(6)));