function plotPath(ks,lambdas,nn,path)
%plot a path in the Brillouin zone
%nn is number of bands
%path should be a mx2 vector, where m is number of points in path

e=1.602e-19;
ts=linspace(0,1,100);
%but vectors in a form interp2 can handle
kxs=unique(ks(:,1))';
kys=unique(ks(:,2));
kxgrid=kron(kxs,ones(length(kys),1));
kygrid=kron(kys,ones(1,length(kxs)));

%inelegant way of solving rounding errors that take pts on edge outside BZ
corr=1-1e-5;

pathx=path(:,1);
pathx(pathx>kxs(end))=kxs(end)*corr;
pathx(pathx<kxs(1))=kxs(1)*corr;

pathy=path(:,2);
pathy(pathy>kys(end))=kys(end)*corr;
pathy(pathy<kys(1))=kys(1)*corr;

path=[pathx,pathy];

energies=[];
for ii=1:nn
    vals=zeros(size(path,1)-1,100);
    n_energies=[];
    lgrid=reshape(lambdas(:,ii),[length(kys),length(kxs)]);
    %define interpolation function for band
    nband=@(kx,ky)interp2(kxgrid,kygrid,lgrid,kx,ky);
    for jj=1:size(path,1)-1
        ptsx=path(jj,1)*(1-ts)+path(jj+1,1)*ts;
        ptsy=path(jj,2)*(1-ts)+path(jj+1,2)*ts;
        vals(jj,:)=nband(ptsx,ptsy);
    end
    
    n_energies=reshape(vals',[1,(size(path,1)-1)*100]);
    energies=cat(1,energies,n_energies);
end
plot(energies'/e,'.');
ylabel('eV');




end

