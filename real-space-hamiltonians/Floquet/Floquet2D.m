%fundamental constants
hbar=1.054e-34; %J.s
m=9.10938291e-31; %kg
c=2.99792458e8; %m/s
e=1.602176565e-19; %C

%Set time periodicity
k_stat_lattice=2*pi/104e-9;
Erecoil=hbar^2*k_stat_lattice^2/(2*m);
omega=1.71*Erecoil;
T=2*pi/omega;
%Define lattice vectors in real space. In m.
a1=-2*pi*[sqrt(3)/2,0.5]/(sqrt(3)*ktime);
a2=2*pi*[sqrt(3)/2,-0.5]/(sqrt(3)*ktime);
amax=max(norm(a1),norm(a2));
unit_cell_area=norm(cross([a1,0],[a2,0]));
%Need a function which is our potential in the unit cell and zero
%outside the unit cell.
xv=[0,a1(1),a1(1)+a2(1),a2(1)];
yv=[0 a1(2), a1(2)+a2(2), a2(2)];
charfn_unit_cell=@(x,y)inpolygon(x,y,xv,yv);
%potential_unit_cell=@(x,y)(potential_form(x,y).*WignerSeitz(x,y,a1,a2));
%Removing offset not useful here.
potential=@(x,y,t)(potential_form(x,y,t)).*charfn_unit_cell(x,y);
%Define potential in Joules
potential_form=@(x,y,t)paper_potential(x,y,0,t); %.*cos(omega*t);  

%Number of Fourier modes, reciprocal lattice vectors, and Brillouin zone q vectors to keep.
%These should all be odd.
Nf=3;
nrecb1=5;
nrecb2=5;
Nbrill1=51;
Nbrill2=51;

%Create fourier modes
nfvals=-(Nf-1)/2:1:(Nf-1)/2;
%Create Reciprocal lattice vectors
[b1,b2]=getRecpBasis(a1,a2);
recb1_n=(-(nrecb1-1)/2:1:(nrecb1-1)/2);
b1RecpVects=kron(transpose(recb1_n),b1);
recb2_n=(-(nrecb2-1)/2:1:(nrecb2-1)/2);
b2RecpVects=kron(transpose(recb2_n),b2);
%create allowed k-vectors in first Brillouin zone.
kb1_n=(-(Nbrill1-1)/2:1:(Nbrill1-1)/2)/Nbrill1;
kb1Vects=kron(transpose(kb1_n),b1);
kb2_n=(-(Nbrill2-1)/2:1:(Nbrill2-1)/2)/Nbrill2;
kb2Vects=kron(transpose(kb2_n),b2);

%define complex exponential function
four_space_exp=@(x,y,k1,k2) exp(1i*(x*k1+y*k2));
four_time_exp=@(t,w) exp(1i*w*t);

%create potential matrix
%matrix indices=dictionary=(-n1,-n2),(-n1,-n2+1),...(-n1,n2),(-n1+1,-n2)...
%this means that the index(l1,l2)=(n1+l1)*(2n2+1)+(l2+n2+1)
%implement this order using kronecker product as below.

recpvects=kron(b1RecpVects,ones(1,length(recb2_n)).')+kron(ones(1,length(recb1_n)).',b2RecpVects);
Umat=zeros(nrecb1*nrecb2,nrecb1*nrecb2,Nf);

for nn=1:Nf
    for ii=1:size(recpvects,1)
        for jj=1:size(recpvects,1)
            Q=recpvects(jj,:)-recpvects(ii,:);
            angfreq=omega*nfvals(nn);
            Umat(ii,jj,nn)=(1/(T*unit_cell_area))*integral3(@(x,y,t)(potential(x,y,t).*four_space_exp(x,y,Q(1),Q(2)).*four_time_exp(t,angfreq)),-amax,amax,-amax,amax,0,T);

        end
    end
end

%create skeleton matrix. This contains terms that do not vary with q.
H_skeleton=kron(diag(nfvals),eye(length(recpvects)))*hbar*omega;
%Blocks for each time fourier component
for nn=1:Nf
    H_skeleton=H_skeleton+kron(diag(ones(1,Nf-abs(nfvals(nn))),nfvals(nn)),Umat(:,:,nn));
end


%create array to hold eigenvalues for each q value.
lambda=[];
kvects=kron(ones(Nbrill2,1),kb2Vects)+kron(kb1Vects,ones(Nbrill1,1));
for ii=1:size(kvects,1)
    k=kvects(ii,:);
    %delta part
    delta_vect=zeros(1,nrecb1*nrecb2);
    for jj=1:size(recpvects,1)
        Q=recpvects(jj,:);
        delta_vect(jj)=(hbar^2/(2*m)).*norm(k-Q)^2;
    end
    delta_part=kron(eye(Nf),diag(delta_vect));
    %Create and diagonalize matrix. Append as next row of lambda
    floquet_H=delta_part+H_skeleton;
    lambda=cat(1,lambda,transpose(eig(floquet_H)));
end

%reduce k's to Brillouin zone
kvalscell=num2cell(kvects,2);
getkBZ=@(pt)reduceToWignerSeitz(pt(1),pt(2),b1,b2);
kvalsBZ=cell2mat(cellfun(getkBZ,kvalscell,'UniformOutput',0));

%plot
scatter3(kvalsBZ(:,1),kvalsBZ(:,2),real(lambda(:,1))/e)
hold on;
scatter3(kvalsBZ(:,1),kvalsBZ(:,2),real(lambda(:,2))/e)
hold on;
scatter3(kvalsBZ(:,1),kvalsBZ(:,2),real(lambda(:,3))/e)
hold on;
scatter3(kvalsBZ(:,1),kvalsBZ(:,2),real(lambda(:,4))/e)
xlabel('Kx')
ylabel('Ky')
zlabel('Energy [eV]')

