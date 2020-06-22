function plotBands(ks,lambdas,jj)
%ks as an nx2 matrix (n is # of points in Brillouin zone)
%lambdas as an nxm matrix (m is the # of bands)
%jj is number of bands to plot

e=1.602e-19;

for ii=1:jj
    scatter3(ks(:,1)*1e-9,ks(:,2)*1e-9,lambdas(:,ii)/e);
    hold on;
end

xlabel('Kx [1/nm]');
ylabel('Ky [1/nm]');
zlabel('Energy [eV]');

end

