%produce 2D plot gamma to M to X to gamma

hbar=1.054e-34; %J.s
m=9.10938291e-31; %kg
c=2.99792458e8; %m/s
e=1.602176565e-19; %C
klatt=0;
Erecoil=hbar^2*klatt^2/(2*m);

%Define lattice vectors in real space, in nm.
a1=1*[1 0];
a2=1*[0 1];
%Find reciprocal lattice basis vectors
[b1,b2]=getRecpBasis(a1,a2);

%Brillouin zone vectors
kdeltaM=cat(2,linspace(0,pi,100)',linspace(0,pi,100)');
kMX=cat(2,pi*ones(100,1),1*linspace(0,pi,100)');
kXdelta=cat(2,linspace(0,pi,100)',zeros(100,1));
kvects=cat(1,kdeltaM,kMX,kXdelta);

tdd=1;
tpp=1;
tcrosspp=0.2;
tpd=1;
delta=1;
sigma=0;


lambdas=[];
eigvects1=[];
eigvects2=[];
eigvects3=[];

Hout=@(kx,ky) [[delta-2*tdd*(cos(kx)+cos(ky)),2*1i*tpd*sin(kx),2*1i*tpd*sin(ky)];[-2*1i*tpd*sin(kx),2*tpp*cos(kx)-2*tcrosspp*cos(ky),1i*sigma];[-2*1i*tpd*sin(ky),-1i*sigma,2*tpp*cos(ky)-2*tcrosspp*cos(kx)]];
%typ in paper. Give H_33=sin(kx). Should be sin(ky).

for ii=1:length(kvects)
    kx=kvects(ii,1);
    ky=kvects(ii,2);
  
    [V,D]=eig(Hout(kx,ky));
    eigvalsvector=[D(1,1),D(2,2),D(3,3)];
    vals=sort(eigvalsvector);
    vects1=V(:,find(eigvalsvector==vals(1)))';
    vects2=V(:,find(eigvalsvector==vals(2)))';
    vects3=V(:,find(eigvalsvector==vals(3)))';
    lambdas=cat(1,lambdas,vals);
    eigvects1=cat(1,eigvects1,vects1);
    eigvects2=cat(1,eigvects2,vects2);
    eigvects3=cat(1,eigvects3,vects3);
end


EdeltaM=real(lambdas(1:100,:));
EMX=real(lambdas(101:200,:));
EXdelta=real(lambdas(201:end,:));

% figure;
% plot(linspace(0,pi,100),EdeltaM,'b.')
% hold on;
% plot(linspace(2*pi,pi,100),EMX,'b.')
% hold on;
% plot(linspace(3*pi,2*pi,100),EXdelta,'b.')
% hold on;
% a=ylim;
% plot([pi,pi],[a(1),a(2)],'black')
% hold on;
% plot([2*pi,2*pi],[a(1),a(2)],'black')
% hold on;
% plot([3*pi,3*pi],[a(1),a(2)],'black')
% 
%Berry curvature
grdHx=@(kx,ky)[[2*tdd*sin(kx), 2*1i*tpd*cos(kx), 0];[-2*1i*tpd*cos(kx),-2*tpp*sin(kx),0];[0,0,2*tcrosspp*sin(kx)]];
grdHy=@(kx,ky)[[2*tdd*sin(ky),0,2*1i*tpd*cos(ky)];[0,2*tcrosspp*sin(ky),0];[-2*1i*tpd*cos(ky),0,-2*tpp*sin(ky)]];
grdHz=@(kx,ky) 0;

kxs=cat(2,linspace(-pi,pi,50)',zeros(50,1));
kys=cat(2,zeros(50,1),linspace(-pi,pi,50)');
kvals=kron(kxs,ones(length(kys),1))+kron(ones(length(kxs),1),kys);
berrycurve=zeros(length(kvals),3);
for ii=1:length(kvals)
    kx=kvals(ii,1);
    ky=kvals(ii,2);
    
    [V,D]=eig(Hout(kx,ky));
    eigvalsvector=[D(1,1),D(2,2),D(3,3)];
    [vals,sort_order]=sort(eigvalsvector);
    vects=V(:,sort_order);
    
    vects(:,1)=vects(:,1)/norm(vects(:,1));
    vects(:,2)=vects(:,2)/norm(vects(:,2));
    vects(:,3)=vects(:,3)/norm(vects(:,3));
    
    
    ipvectone=zeros(3,1);
    ipvecttwo=zeros(3,1);
    pertsum=0;
    
    nn=2; %band to compute curvature for
    bands=[1,2,3];
    for jj=bands(bands~=nn)
        ipvectone(1)=vects(:,nn)'*grdHx(kx,ky)*vects(:,jj);
        ipvectone(2)=vects(:,nn)'*grdHy(kx,ky)*vects(:,jj);
        ipvectone(3)=0;
        ipvectone=ipvectone/(eigvalsvector(nn)-eigvalsvector(jj));
        
        ipvecttwo(1)=vects(:,jj)'*grdHx(kx,ky)*vects(:,nn);
        ipvecttwo(2)=vects(:,jj)'*grdHy(kx,ky)*vects(:,nn);
        ipvecttwo(3)=0;
        ipvecttwo=ipvecttwo/(eigvalsvector(nn)-eigvalsvector(jj));
        
        pertsum=pertsum+cross(ipvectone,ipvecttwo);
    end
    berrycurve(ii,:)=-imag(pertsum)';
end

figure;
quiver(kvals(:,1),kvals(:,2),berrycurve(:,1),berrycurve(:,2))